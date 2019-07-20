#include <iostream>
#include <ctime>    
#include <cstdlib>  
//#include "C:\\Users\\\baron\\Desktop\\School\\Research\\Fa 2018 - Sum 2019\\(1)Code\\C++ Code Directory\\PoscarInfo.h"
//#include "C:\\Users\\\baron\\Desktop\\School\\Research\\Fa 2018 - Sum 2019\\(1)Code\\C++ Code Directory\\Calculations.h"
//#include "C:\\Users\\\baron\\Desktop\\School\\Research\\Fa 2018 - Sum 2019\\(1)Code\\C++ Code Directory\\AtomSwitcher.h"
//#include "C:\\Users\\\baron\\Desktop\\School\\Research\\Fa 2018 - Sum 2019\\(1)Code\\C++ Code Directory\\NPA.h"
#include "PoscarInfo.h"
#include "Calculations.h"
//#include "AtomSwitcher.h"
#include "NPA.h"

std::string infileLoc = DEFAULT_POSCAR_OPEN_PATH;
std::string bondInfoLoc = DEFAULT_BONDINFO_OPEN_PATH;
std::string writeLoc = (DEFAULT_POSCAR_OPEN_PATH + "new").c_str();
std::string outputDataLoc = ((PATH_TO_FILES + "AutoAtomSwitcherInfo.csv").c_str());
int nAttempts = 50; ///number of attempts to retry making a defected POSCAR before code gives up and moves to the next one


//STRUCTS, FUNCTION DECL, ETC---------------------------------------------------------------------------------------------------
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

struct idToSwitch
{
	int id_;
	int count_; ///num of this specific ID
	std::string switchTo;
	bool okToSwitch; ///for stopping specific atoms from being switched - for example, maybe user wants to change all SiC to SiSi except cases where a Si is bonded to H
	bool copied;

	idToSwitch(int id__, std::string s)
	{
		id_ = id__;
		switchTo = s;
		count_ = 1;
		copied = false;
		okToSwitch = true;
	}
};

struct notOkIDs
{
	int id_;
	bool ok;

	notOkIDs(int _id_)
	{
		id_ = _id_;
		ok = false;
	}
};

std::vector <bondInfo> bondInfoVect;
bool goodBondInfo = false;

struct moveInfo
{
	Coords notHyd, hyd;
	double necessaryDisp;
	long double unitVectHyd [3];

	moveInfo()
	{
		notHyd;
		hyd;
		necessaryDisp = -99;
	}

	moveInfo(Coords hyd_, Coords notHyd_, std::string swapNotHydTo)
	{
		notHyd = notHyd_;
		hyd = hyd_;

		notHyd.atomType = swapNotHydTo;

		//Find necessary displacement
		if (bondInfoVect.size() == 0)
			std::cout << "constructor moveInfo(): bondInfoVect size is 0";
		for (int i = 0; i < bondInfoVect.size(); i++)
		{
			goodBondInfo = false;
			if ((new atomPair(notHyd, hyd))->type == bondInfoVect[i].type1 || (new atomPair(notHyd, hyd))->type == bondInfoVect[i].type2)
			{
				necessaryDisp = bondInfoVect[i].physBondDist;
				goodBondInfo = true;
				break;
			}
		}

		//Get the unit vector for hydrogen
		long double displacements[3] = { hyd.a - notHyd.a, hyd.b - notHyd.b, hyd.c - notHyd.c}; ///displacements from hydrogen to not hydrogen
		for (int i = 0; i < 3; i++)
			unitVectHyd[i] = displacements[i] / dist_(hyd, notHyd);
	}
};

struct idsAndCounts
{
	bool hyd;
	bool considered = false;
	int id;
	int Bondcount;

	idsAndCounts(Coords crds)
	{
		if (crds.atomType == "H")
			hyd = true;
		else
			hyd = false;

		id = crds.id;
	}
};

int main()
{
	srand(time(NULL)); ///seed RNG

	//READ IN AND CHECK INFILE POSCAR-----------------------------------------------------------------------------------------------
			//-------------------------------------------------------------------------------------------------------------------------------
	Poscar P("readAll", infileLoc.c_str());
	Poscar copyFormat = P;
	P.fetchAtomBonds(bondInfoLoc.c_str());
	P.atomBonds.pop_back();
	copyFormat.atomBonds = P.atomBonds;
	//Checks for missing atom IDs in bonds (each atom should have at least one bond)
	{
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
				return 0;
			}
		}
	}
	//Checks for non-physical species
	{
		std::cout << "Input file: ";
		bool err = false;

		npaPoscar W("readAll", infileLoc.c_str());
		W.atomBonds = P.atomBonds;
		W.atomCoords = P.atomCoords;
		W.fetchThisSpeciesVector(bondInfoLoc.c_str());
		for (int i = 0; i < W.thisSpeciesVect.size(); i++)
		{
			if (W.thisSpeciesVect[i].mainAtom.atomType == "H")
				if (W.thisSpeciesVect[i].atomBonds_.size() != 1)
				{
					std::cout << "\nError: Hydrogen " << W.thisSpeciesVect[i].mainAtom.id + 1 << " has " << W.thisSpeciesVect[i].atomBonds_.size() << " bonds, but it should have 1 bond";
					err = true;
				}
			if (W.thisSpeciesVect[i].mainAtom.atomType != "H")
				if (W.thisSpeciesVect[i].atomBonds_.size() != 4)
				{
					std::cout << "\nError: " << W.thisSpeciesVect[i].mainAtom.atomType << " " << W.thisSpeciesVect[i].mainAtom.id + 1 << " has " << W.thisSpeciesVect[i].atomBonds_.size() << " bonds, but it should have 4 bonds";
					err = true;
				}
		}

		if (err)
		{
			std::cout << "\n";
			return 0;
		}
		else
			std::cout << "ok\n";
	}
	//Checks bond info file for anything other than 6 entries (this is very specific but this code is only for SiC anyways)
	{
		if (P.atomTypes.size() == 3)
		{
			std::ifstream infile(bondInfoLoc.c_str());
			int sum = -1;
			std::string garbage;
			while (!infile.eof())
			{
				std::getline(infile, garbage);
				sum++;
			}

			if (sum != 6)
			{
				std::cout << "badly formatted bond data file (nLines = " << sum << ".  Was expecting 6.\n";
				return 0;
			}
			infile.close();
		}
	}

	//---------------------------------------------------------------------------------------------------------------------------------

	//OUTPUT INFO ABOUT INFILE, GRAB USER INPUT-------------------------------------------------------------------------------------
			//Bond information
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

	//User input
	//Prompt for what kind of change to make
	char choice;
	std::string tmpType1, tmpType2; ///switch from
	std::string tmpType1_, tmpType2_; ///switch to
	int loA, hiA, loB, hiB, nModels;
	std::cout << "----------\n(a) Switch C-H to Si-H\n(b) Switch Si-H to C-H\n(c) Number of models to randomly generate\n";
	std::cout << "Enter values in the following format:\n<a low> <a high>\n<b low> <b high>\n<c>\n";
	std::cin >> loA >> hiA >> loB >> hiB >> nModels;

	//Flip low and high nums in case user dosent know how to read
	if (loA > hiA)
	{
		int tmp = loA;
		loA = hiA;
		hiA = tmp;
	}
	if (loB > hiB)
	{
		int tmp = loB;
		loB = hiB;
		hiB = tmp;
	}

	//Generates outfile, write header-----------------------------------------------------------------------------------------------
				std::ofstream codeData(outputDataLoc.c_str());
	codeData << "Lower Bound C-H to Si-H,Upper Bound C-H to Si-H,Lower Bound Si-H to C-H,Upper Bound Si-H to C-H,Max number of attempts per file\n"
		<< loA << "," << hiA << "," << loB << "," << hiB << "," << nAttempts << "\n";
	codeData << "Name,Defected? (0 = n),Number of attempts to make,Num C-H to Si-H,Num Si-H to C-H\n";

	//MAIN LOOP:  ONE ITERATION FOR EACH SUCESSFUL POSCAR FILE GENERATED------------------------------------------------------------
	for (int mainLoopNum = 0; mainLoopNum < nModels; mainLoopNum++)
	{
		//Get the number of switches to make for each type for a certian output----------------------------------------------------
		int aSwitches = (rand() % (hiA - loA + 1)) + loA; ///CH to SiH
		int bSwitches = (rand() % (hiB - loB + 1)) + loB; ///SiH to CH

		int thisAttemptNum;
		bool defected = false;

		//Defect loop: Continues to loop until a non-defect POSCAR is made, or the number of allowed attempts runs out ------------
		for (thisAttemptNum = 0;; thisAttemptNum++)
		{
			//Exit conditions
			if (thisAttemptNum >= nAttempts)
			{
				defected = true;
				break;
			}

			//Defect info initialization
			std::vector <std::vector <Coords> > defInfoVect;

			//Create new instance of POSCAR initially copying the input file---------------------------------------------------------
			Poscar R = P;
			//Switch loop: two iterations, one for each type of switch to be made-----------------------------------------------------
			for (int switchType = 0; switchType < 2; switchType++)
			{
				//Choosing what kind of type to switch: 0 corresponds to type A, 1 to type B.  Choosing num to switch------------------
				if (switchType == 0)
					choice = 'a';
				if (switchType == 1)
					choice = 'b';

				int numToSwitch = 0;
				if (choice == 'a')
				{
					tmpType1 = "CH";
					tmpType2 = "HC";
					tmpType1_ = "SiH";
					tmpType2_ = "HSi";
					numToSwitch = aSwitches;
					
				}
				else if (choice == 'b')
				{
					tmpType1 = "SiH";
					tmpType2 = "HSi";
					tmpType1_ = "CH";
					tmpType2_ = "HC";
					numToSwitch = bSwitches;
				}
				//Seperating the types into two seperate strings to ID atoms by----------------------------------------------------------
								std::string type1 = "", type2 = "";
				bool one = true;
				std::string  _tmpType1_ = (tmpType1_ + "!").c_str();
				for (int i = 0; i < _tmpType1_.size(); i++)
				{
					if (one)
						type1 = (type1 + _tmpType1_[i]).c_str();
					else
						type2 = (type2 + _tmpType1_[i]).c_str();

					if (isupper(_tmpType1_[i + 1]))
						one = false;
					if (_tmpType1_[i + 1] == '!')
						break;
				}
				//Mark atoms that have already been switched in a pervious loop, so as to not switch them again.  Initialize list of atoms that are OK to switch---------------------
				std::vector <notOkIDs> notOk;
				///Cases where an atom has already been counted
				for (int i = 0; i < R.atomCoords.size(); i++)
					if (R.atomCoords[i].extraVal == true) ///if an atom has been copied
						notOk.push_back(*new notOkIDs(R.atomCoords[i].id));
				//Ititialize bond list with IDs that are ok to switch
				for (int i = 0; i < R.atomBonds.size(); i++)
				{
					for (int n = 0; n < 2; n++)
						for (int j = 0; j < notOk.size(); j++)
							if (R.atomBonds[i].pairedAtoms[n].id == notOk[j].id_)
								R.atomBonds[i].pairedAtoms[n].extraInfo = "no";
				}
				//Mark atoms to switch, put them in vector of switches to make---------------------------------------------------------
								std::vector <idToSwitch> switches;
				for (int i = 0; i < R.atomBonds.size(); i++)
					if (R.atomBonds[i].type == tmpType1 || R.atomBonds[i].type == tmpType2)
					{
						//Specific special cases
						if (choice == 'a')
							for (int n = 0; n < 2; n++)
								if (R.atomBonds[i].pairedAtoms[n].atomType == "C" && R.atomBonds[i].pairedAtoms[n].extraInfo != "no")
									switches.push_back(*new idToSwitch(R.atomBonds[i].pairedAtoms[n].id, "Si"));

						if (choice == 'b')
							for (int n = 0; n < 2; n++)
								if (R.atomBonds[i].pairedAtoms[n].atomType == "Si" && R.atomBonds[i].pairedAtoms[n].extraInfo != "no")
									switches.push_back(*new idToSwitch(R.atomBonds[i].pairedAtoms[n].id, "C"));
					}
				if (switches.size() == 0)
					break;
				//Get rid of duplicate switches / initialize counts----------------------------------------------------------------------
								std::vector <idToSwitch> tmp;
				for (int i = 0; i < switches.size(); i++)
				{
					for (int j = 0; j < switches.size(); j++)
						if (switches[i].id_ == switches[j].id_ && i != j && !switches[i].copied)
						{
							switches[i].count_++;
							switches[j].copied = true;
						}
					if (!switches[i].copied)
						tmp.push_back(switches[i]);
				}
				switches = tmp;
				//Make switches for this switch type--------------------------------------------------------------------------------------
								std::random_shuffle(switches.begin(), switches.end());
				std::random_shuffle(R.atomCoords.begin(), R.atomCoords.end());
				for (int i = 0, breakCond = 0; i < R.atomCoords.size() && breakCond < numToSwitch; i++)
					for (int j = 0; j < switches.size(); j++)
					{
						if (R.atomCoords[i].id == switches[j].id_)
						{
							R.atomCoords[i].atomType = switches[j].switchTo;
							R.atomCoords[i].extraVal = true; ///flag to show that the atom has been considered already
							breakCond += switches[j].count_; ///so switching a single atom with n bonds counts as switching n bonds, instead of just one
						}
					}
				//Update Poscar-----------------------------------------------------------------------------------------------------------
								R.removeDuplicates();
				R.updateAll();
			}

			//Move hydrogen atoms closer or further away from switched atoms---------------------------------------------------------------
						R.convertToCartesian();

			//Fill bond info vector
			std::ifstream infile(bondInfoLoc.c_str());
			while (!infile.eof())
			{
				std::string tmp1, tmp2;
				double dist;
				infile >> tmp1 >> tmp2 >> dist;

				bondInfoVect.push_back(*new bondInfo(tmp1, tmp2, dist)); 
			}

			//Check the original bonds that had hydrogens in them
			for (int i = 0; i < R.atomBonds.size(); i++)
			{
				int idToCompare = -99;
				int thisHid = -99;
				Coords notHydrogen;
				if (R.atomBonds[i].pairedAtoms[0].atomType == "H")
				{
					thisHid = R.atomBonds[i].pairedAtoms[0].id;
					idToCompare = R.atomBonds[i].pairedAtoms[1].id;
					notHydrogen = R.atomBonds[i].pairedAtoms[1];
				}
				else
					if (R.atomBonds[i].pairedAtoms[1].atomType == "H")
					{
						thisHid = R.atomBonds[i].pairedAtoms[1].id;
						idToCompare = R.atomBonds[i].pairedAtoms[0].id;
						notHydrogen = R.atomBonds[i].pairedAtoms[0];
					}

				///If an instance of a hydrogen bond is found
				if (idToCompare != -99)
					//go through atomCoords to make sure it has been swapped (extraVal == true)
					for (int j = 0; j < R.atomCoords.size(); j++)
						if (R.atomCoords[j].extraVal == true && idToCompare == R.atomCoords[j].id) ///if it has been swapped...
						{
							//Edit the hydrogen's coords:
							///first set up the move info (an odd way of doing this, but I plan on using the moveInfo struct in later places anyways)
							moveInfo move;
							std::string swapTo;
							if (R.atomBonds[i].pairedAtoms[0].atomType == "H")
							{
								if (R.atomBonds[i].pairedAtoms[1].atomType == "C")
									swapTo = "Si";
								else
									swapTo = "C";
								move = *new moveInfo(R.atomBonds[i].pairedAtoms[0], R.atomBonds[i].pairedAtoms[1], swapTo);
							}
							else
							{
								if (R.atomBonds[i].pairedAtoms[0].atomType == "C")
									swapTo = "Si";
								else
									swapTo = "C";
								move = *new moveInfo(R.atomBonds[i].pairedAtoms[1], R.atomBonds[i].pairedAtoms[0], swapTo);
							}

							///then change the hydrogen's coordinates
							for (int k = 0; k < R.atomCoords.size(); k++)
							{										//vvv just to make absolutly sure im changing a hydrogen
								if (R.atomCoords[k].id == thisHid && R.atomCoords[k].atomType == "H") //once the correct hydrogen is found
								{
									R.atomCoords[k].a = notHydrogen.a + (move.unitVectHyd[0] * (move.necessaryDisp - 0.01));
									R.atomCoords[k].b = notHydrogen.b + (move.unitVectHyd[1] * (move.necessaryDisp - 0.01));
									R.atomCoords[k].c = notHydrogen.c + (move.unitVectHyd[2] * (move.necessaryDisp - 0.01));
								}
							}

						}
			}
			//Copy format, write new file (which possibly has defects)---------------------------------------------------------------------------
						R.copyFormatting(copyFormat);
			R.write((writeLoc + std::to_string(mainLoopNum)).c_str());
			//Check to make sure newly created Poscar has no defects.  If ok, break from defect loop. Otherwise reset the defect loop---------
			bool err = false;
			{
				npaPoscar W("readAll", (writeLoc + std::to_string(mainLoopNum)).c_str());
				W.fetchAtomBonds(bondInfoLoc.c_str());
				W.atomBonds.pop_back(); ///this wasnt in the original code (that worked, mind you) but I think it should have been.  If something fails, try removing this
				W.fetchThisSpeciesVector(bondInfoLoc.c_str());
				for (int i = 0; i < W.thisSpeciesVect.size(); i++)
				{
					if (W.thisSpeciesVect[i].mainAtom.atomType == "H")
						if (W.thisSpeciesVect[i].atomBonds_.size() != 1)	
						{
							err = true;
							std::vector <Coords> defInfo;
							defInfo.push_back(W.thisSpeciesVect[i].mainAtom);
							for (int j = 0; j < W.thisSpeciesVect[i].atomBonds_.size(); j++)
								defInfo.push_back(W.thisSpeciesVect[i].atomBonds_[j].pairedAtoms[1]); ///paired atoms [1] is always not the central atom
							defInfoVect.push_back(defInfo);
						}
					if (W.thisSpeciesVect[i].mainAtom.atomType != "H")
						if (W.thisSpeciesVect[i].atomBonds_.size() != 4)
						{
							err = true;
							std::vector <Coords> defInfo;
							defInfo.push_back(W.thisSpeciesVect[i].mainAtom);
							for (int j = 0; j < W.thisSpeciesVect[i].atomBonds_.size(); j++)
								defInfo.push_back(W.thisSpeciesVect[i].atomBonds_[j].pairedAtoms[1]); ///paired atoms [1] is always not the central atom
							defInfoVect.push_back(defInfo);
						}
				}
			}
			if (!err)
			{
				std::ofstream defectOut_;
				defectOut_.open(("POSCAR" + std::to_string(mainLoopNum) + "defectInfo").c_str(), std::ios_base::app);
				defectOut_ << "Completed\n";
				break;
			}
			//append defect info to this poscars defect text file
						else
			{
				std::ofstream defectOut;
				defectOut.open(("POSCAR" + std::to_string(mainLoopNum) + "defectInfo").c_str(), std::ios_base::app);
				defectOut << "\nAttempt ( " << thisAttemptNum + 1 << " ): [Total Defects: " << defInfoVect.size() << "]----------------------\n";
				for (int i = 0; i < defInfoVect.size(); i++)
				{
					defectOut << "Main Atom: " << defInfoVect[i][0].atomType << defInfoVect[i][0].id << "\n";
					defectOut << "Coords: < " << defInfoVect[i][0].a << " , " << defInfoVect[i][0].b << " , " << defInfoVect[i][0].c << " >\n";
					for (int j = 1; j < defInfoVect[i].size(); j++)
					{
						defectOut << "\tBonded Atom ( " << j << " ): " << defInfoVect[i][j].atomType << defInfoVect[i][j].id << " || Dist from Central: " << dist_(defInfoVect[i][0], defInfoVect[i][j]) << "\n";
						defectOut << "\t-Coords: < " << defInfoVect[i][j].a << " , " << defInfoVect[i][j].b << " , " << defInfoVect[i][j].c << " >\n";
					}
				}
			}
		}

		std::cout << "File " << mainLoopNum;
		if (!defected)
			std::cout << " completed successfully\n";
		else
			std::cout << " could not be written without defects in the given number of attempts\n";
		//Write generation information to csv----------------------------------------------------------------------------------------------
		codeData << (writeLoc + std::to_string(mainLoopNum)).c_str() << "," << defected << "," << thisAttemptNum + 1 << "," << aSwitches << "," << bSwitches << "\n";
	}

	codeData.close();
	return 0;
}
