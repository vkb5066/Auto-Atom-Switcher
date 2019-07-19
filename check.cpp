#include <iostream>
#include "PoscarInfo.h"
#include "Calculations.h"

std::string infileLoc = DEFAULT_POSCAR_OPEN_PATH;
std::string bondInfoLoc = DEFAULT_BONDINFO_OPEN_PATH;
std::string writeLoc = (DEFAULT_POSCAR_OPEN_PATH + "new").c_str();

struct bondSwitchAssist
{
        int code;
        int count;
        std::string type1, type2;
        std::vector <atomPair> atomBonds_;

        bondSwitchAssist(int cd, std::vector <atomPair> atmprVect)
        {
                code = cd;
                count = atmprVect.size();
                type1 = atmprVect[0].type;
                type2 = atmprVect[0].type2;
                atomBonds_ = atmprVect;
        }
};

int main()
{
	Poscar P("readAll", (infileLoc + "new").c_str());
	P.fetchAtomBonds(bondInfoLoc.c_str());
	P.atomBonds.pop_back();

	//Search for missing IDs from the list of bonds compared to the atoms being read in
	for (int i = 0; i < P.atomCoords.size(); i++)
	{
		bool missingAtomId = true;
		for (int j = 0; j < P.atomBonds.size(); j++)
			for (int n = 0; n < 2; n++)
				if (P.atomCoords[i].id == P.atomBonds[j].pairedAtoms[n].id)
					missingAtomId = false;
		if (missingAtomId)
		{
			std::cout << "At least one atom from the input file is not bonded to anything, by the definitions given in the bond data info file\n";
			//system ("pause");
			return 0;
		}
	}

	//Output current bond Information
	std::vector <bondSwitchAssist> bondSwitch;
	std::cout << "Input file information:\n";
	for (int i = 0; i < P.atomBonds.size(); i++)
	{

		if (P.atomBonds[i].extraInfo != "counted")
		{
			int count = 0;
			static int k = 0;
			std::vector <atomPair> tmp;
			for (int j = 0; j < P.atomBonds.size(); j++)
			{
				if ((P.atomBonds[i] ^= P.atomBonds[j]) && P.atomBonds[j].extraInfo != "counted")
				{
					P.atomBonds[j].extraInfo = "counted";
					count++;
					tmp.push_back(P.atomBonds[j]);
				}
			}
			k++;
			std::cout << "(" << k << ") " << P.atomBonds[i].type << " : " << count << "\n";
			bondSwitch.push_back(*new bondSwitchAssist(k, tmp));
		}

	}

	return 0;
}
