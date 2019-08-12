#!/bin/csh -f

if (-d poscars/) then
	rm -r poscars
endif

mkdir poscars
mkdir poscars/defectInfo

#RUN EXAMPLE:
#in order: write progress to screen, create defect log, create csv with info about POSCARS, be lenient with what is considered a defect
#./drive verbose log csv len

./drive verbose log csv


#moving stuff around to declutter run directory
if (-f verb) then
	mv -v *new* poscars/.
	mv -v *defect* poscars/defectInfo/.
	rm verb
else
	mv *new* poscars/.
	mv *defect* poscars/defectInfo/.
endif

echo "Complete"
