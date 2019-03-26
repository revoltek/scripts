import casacore.tables as tab
a=tab.table('/export/scratch/AG_deGasperin/alex/VLA/PSZ1_G096.89+24.17/fdg/lustre/aoc/ftp/e2earchive/stage/BI0004/obs_3_test.MS',readonly=False)
a.unlock()
tc=a.col('DATA')
print(tc[14773200-1][0])