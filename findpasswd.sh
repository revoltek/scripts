#!/bin/bash
# findpasswd.sh "search pattern"
# to add a passwd: gpg file.gpg
# and then: gpg -c file

gpg ~/altro/memo.gpg
if [ $? -ne 0 ]; then
	echo "ERROR: wrong code?"
	exit 1
fi
if [ ${1}x == ''x ]; then
	cat ~/altro/memo
else
	grep -i $1 ~/altro/memo
fi
rm ~/altro/memo
