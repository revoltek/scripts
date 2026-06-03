# use:
# define: bias, contrast, clip_min, clip_max and scale_type
#
#     my_stretch = CartaFullShaderStretch(bias, contrast, scale_type='SQRT')
#     norm = ImageNormalize(
#        data, 
#        vmin=clip_min, 
#        vmax=clip_max, 
#        stretch=my_stretch, 
#        clip=True
#    )
#    im = ax.imshow(clean_data, origin="lower", cmap="afmhot", norm=norm)


import numpy as np
from astropy.visualization import BaseStretch

class CartaFullShaderStretch(BaseStretch):
    def __init__(self, bias, contrast, scale_type='SQRT', alpha=1000.0, gamma=2.2):
        super().__init__()
        self.uBias = bias
        self.uContrast = contrast
        self.uScaleType = scale_type.upper()
        self.uAlpha = alpha
        self.uGamma = gamma

    def error_function(self, x, c, x0):
        # Implementation of your provided logistic function:
        # y = exp(c * (x - x0)) / (exp(c * (x - x0)) + 1.0)
        # Note: We use np.exp for array compatibility
        y = np.exp(c * (x - x0))
        return y / (y + 1.0)

    def __call__(self, values, clip=True, out=None):
        x = np.array(values, copy=True)

        # --- PHASE 1: Scale Type (BEFORE Bias/Contrast) ---
        if self.uScaleType == 'SQUARE':
            x = x * x
        elif self.uScaleType == 'SQRT':
            x = np.sqrt(x)
        elif self.uScaleType == 'LOG':
            x = np.log(self.uAlpha * x + 1.0) / np.log(self.uAlpha + 1.0)
        elif self.uScaleType == 'POWER':
            x = (np.power(self.uAlpha, x) - 1.0) / (self.uAlpha - 1.0)
        elif self.uScaleType == 'GAMMA':
            x = np.power(x, self.uGamma)
        # Default is Linear (no change to x)

        if self.uContrast <= 1.0:
            # --- Linear Path ---
            smoothedBias = 0.5 - self.uBias / 2.0
            x = (x - smoothedBias) * self.uContrast + smoothedBias
        else:
            # --- Sigmoidal Path (Logistic) ---
            smoothedBias = self.uBias / 2.0 + 0.5
            smoothedContrast = self.uContrast - 1.0
            # [1, 2] map to [0, 1] then scaled by 12.0
            smoothedContrast = 0.001 if smoothedContrast == 0.0 else smoothedContrast * 12.0
            
            offset = self.error_function(0.0, smoothedContrast, smoothedBias)
            denominator = self.error_function(1.0, smoothedContrast, smoothedBias) - offset
            
            if denominator <= 0.0:
                denominator = 0.1
                
            x = (self.error_function(x, smoothedContrast, smoothedBias) - offset) / denominator
            
        # Final clamping and Square Root scaling
        x = np.clip(x, 0.0, 1.0)
        
        if out is not None:
            out[:] = x
            return out
        return x

    @property
    def inverse(self):
        return CartaFullShaderInverse(self.uBias, self.uContrast)

class CartaFullShaderInverse(BaseStretch):
    def __init__(self, bias, contrast, scale_type='SQRT', alpha=1000.0, gamma=2.2):
        super().__init__()
        self.uBias = bias
        self.uContrast = contrast
        self.uScaleType = scale_type.upper()
        self.uAlpha = alpha
        self.uGamma = gamma

    def logistic(self, x, c, x0):
        # Forward logistic for re-calculating offset/denominator
        arg = c * (x - x0)
        y = np.exp(np.clip(arg, -100, 100))
        return y / (y + 1.0)

    def __call__(self, values, clip=True, out=None):
        # We start with 'values' which are the final pixel colors [0, 1]
        x = np.array(values, copy=True)

        # --- PHASE 1: Undo Smoothing (Bias/Contrast) ---
        if self.uContrast <= 1.0:
            # Reverse: x = (x - smoothedBias) * uContrast + smoothedBias
            smoothedBias = 0.5 - self.uBias / 2.0
            # Guard against division by zero if contrast is 0
            c_factor = self.uContrast if self.uContrast != 0 else 0.001
            x = ((x - smoothedBias) / c_factor) + smoothedBias
        else:
            # Reverse: Sigmoidal Path
            smoothedBias = self.uBias / 2.0 + 0.5
            smoothedContrast = (self.uContrast - 1.0) * 12.0
            if smoothedContrast <= 0: smoothedContrast = 0.001

            offset = self.logistic(0.0, smoothedContrast, smoothedBias)
            denominator = self.logistic(1.0, smoothedContrast, smoothedBias) - offset
            
            # Undo normalization: y = (logistic(x) - offset) / denominator
            y_targeted = x * denominator + offset
            # Clip to avoid log(0) or log(inf) in the logit function
            y_targeted = np.clip(y_targeted, 1e-7, 1.0 - 1e-7)
            
            # Logit function (Inverse of Logistic)
            x = (np.log(y_targeted / (1.0 - y_targeted)) / smoothedContrast) + smoothedBias

        # --- PHASE 2: Undo Scale Type ---
        if self.uScaleType == 'SQUARE':
            # Inverse of x^2 is sqrt(x)
            x = np.sqrt(np.clip(x, 0, None))
        elif self.uScaleType == 'SQRT':
            # Inverse of sqrt(x) is x^2
            x = x * x
        elif self.uScaleType == 'LOG':
            # Inverse of log(ax + 1) / log(a + 1)
            # x_orig = (exp(x * log(a + 1)) - 1) / a
            x = (np.exp(x * np.log(self.uAlpha + 1.0)) - 1.0) / self.uAlpha
        elif self.uScaleType == 'POWER':
            # Inverse of (a^x - 1) / (a - 1)
            # x_orig = log(x * (a - 1) + 1) / log(a)
            x = np.log(x * (self.uAlpha - 1.0) + 1.0) / np.log(self.uAlpha)
        elif self.uScaleType == 'GAMMA':
            # Inverse of x^gamma is x^(1/gamma)
            x = np.power(np.clip(x, 0, None), 1.0 / self.uGamma)

        if clip:
            return np.clip(x, 0.0, 1.0, out=out)
        return x
