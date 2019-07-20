#!/bin/csh -f

if (-d poscars/) then
	rm -r poscars
endif

mkdir poscars
mkdir poscars/defectInfo

./drive

mv -v *new* poscars/.
mv -v *defect* poscars/defectInfo/.

echo "complete"
