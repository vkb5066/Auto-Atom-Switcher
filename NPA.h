#pragma once
#include "PoscarInfo.h"
#include "Calculations.h"

//Globals
int MAXBONDS = 4;

//Comparison rule for sorting strings
bool npaStringSort(std::string a, std::string b)
{
	return a < b;
}

//slightly extended version of Coords from the PoscarInfo header
struct npaCoords : public Coords
{
	int count;

	//default constructor
	npaCoords()
	{
		Coords();
		count = -1;
	}

	//Operators
	npaCoords& operator= (const npaCoords &npacoords);
	bool operator== (const npaCoords &) const;
	bool operator!= (const npaCoords &) const;
};
npaCoords& npaCoords::operator= (const npaCoords &npacoords)
{
	Coords::operator=(npacoords);
	//extra stuff here
	return *this;
}
bool npaCoords::operator== (const npaCoords & npacoords) const
{
	return Coords::operator==(npacoords); ///dont need any extra stuff here
}
bool npaCoords::operator!= (const npaCoords &npacoords) const
{
	return Coords::operator!=(npacoords);
}

//struct for below struct's function 'writeString' to make life easier; probably a way to use 'writestring()' without this, but it would be annoying
struct atomTypesAndCounts
{
	std::string type;
	int count;

	///default contructor
	atomTypesAndCounts()
	{
		type = "undf : aotmTypesAndCounts - Default constructor";
		count = -1;
	}

	///constructor that initializes 'type' with the string passed to it
	atomTypesAndCounts(std::string _str_)
	{
		type = _str_;
		count = 0;
	}

	///constructor for type and count
	atomTypesAndCounts(int num, std::string _str_)
	{
		count = num;
		type = _str_;
	}
};

//sort for below class
bool npaBondsSortByType(atomTypesAndCounts a, atomTypesAndCounts b)
{
	return a.type < b.type;
}

//Species struct, main variable type for this class.
///TODO: write operators for the species class
struct Species
{
	Coords mainAtom; ///Coord info about the main (central) atom of the species -- the one every other atom connects to
	std::vector <atomPair> atomBonds_; ///bond info for every bond the main atom is involved in
	std::vector <atomTypesAndCounts> bondedTo;
	int count;
	
	std::string extraInfo; ///just in case

	//default constructor; not very useful, but may be necessary
	Species()
	{
		mainAtom; ///uses default constructor from coords
		count = -1;
		///does not call atomPair constructor because it is possible that this 'Species' is not bonded to anything
		extraInfo = "nada";
	}

	//string constructor; meant to provide some info about the species from reading in the 'speciesStr' string
	Species(std::string speciesString)
	{
		std::vector <int> tmpAtomCounts;
		std::vector <std::string> tmpAtomTypes;

		///Get main atom
		mainAtom.atomType = speciesString.substr(0, speciesString.find(':'));

		for (int i = 0; i < speciesString.size(); i++)
		{
			///End condition
			if (speciesString[i] == '*')
				break;
			///Get info about the rest of the atoms
			else
			{
				//If the current char is a >, then the previous char was an atom count
				if (speciesString[i] == '>')
					tmpAtomCounts.push_back(speciesString[i - 1] - '0'); ///WARNING:  assumes that the atom nmber (number of bonds) is < 10.  Which seems reasonable, but you never know
				//If current char is a !, then the pevious (one or two) char(s) were the atom type
				if (speciesString[i] == '!')
				{
					if (isalpha(speciesString[i - 2])) ///if two spaces back is a letter
					{
						std::string tmp;
						tmp += speciesString[i - 2];
						tmp += speciesString[i - 1];		///WARNING:  this relies heavily on the assumption that nBonds < 10
						tmpAtomTypes.push_back(tmp);
					}
					else
						tmpAtomTypes.push_back(*new std::string = speciesString[i - 1]); ///else one back was the letter
				}
			}
		}

		//Set up some info for the atomBonds_ vector; forcing each atomPair's pairedAtoms[0] position to be the central atom
		for (int i = 0; i < tmpAtomCounts.size(); i++) ///atomCounts and AtomTypes should be equal, so it dosent matter which one is used for max index
		{
			npaCoords tmp;
			tmp.atomType = tmpAtomTypes[i];
			tmp.count = tmpAtomCounts[i];
			atomBonds_.push_back(*new atomPair(mainAtom, tmp));
			bondedTo.push_back(*new atomTypesAndCounts(tmpAtomCounts[i], tmpAtomTypes[i]));
		}

	}

	//Coords constructor; creates a Species with the Coords type provided as the mainAtom
	Species(Coords crds)
	{
		mainAtom = crds;
		extraInfo = "main initialized";
	}

	//function to write the current Species to a string
	std::string writeString(std::vector <std::string> atomTypes)
	{
		std::string tmp;
		std::vector <atomTypesAndCounts> mainBondedTo;

		///need to do this to make sure the order of atoms is correct
		std::vector <std::string> tmpVect = atomTypes;					//TODO:  i've used this alphabetical sorting thing 3 times so far, may be worth
		std::sort(tmpVect.begin(), tmpVect.end(), npaStringSort);		//adjusting PoscarInfo.h to automatically sort atomTypes in the constructor
																		//so this isnt necessary.  Would need to do it carefully, though
		///setting the initial string of potential bonds, with all types in place but all counts set to zero
		for (int i = 0; i < tmpVect.size(); i++)
			mainBondedTo.push_back(*new atomTypesAndCounts(tmpVect[i]));

		tmp = (mainAtom.atomType + ":").c_str();

		///setting counts for each type
		for (int i = 0; i < atomBonds_.size(); i++)
			for (int j = 0; j < mainBondedTo.size(); j++)
				for (int n = 0; n < 2; n++)
				if (atomBonds_[i].pairedAtoms[n].atomType == mainBondedTo[j].type && atomBonds_[i].pairedAtoms[n] != mainAtom)
					mainBondedTo[j].count++;

		///actually create the string
		for (int i = 0; i < mainBondedTo.size(); i++)
		{
			std::stringstream tmpCountStr;
			tmpCountStr << mainBondedTo[i].count;
			tmp = (tmp + tmpCountStr.str() + ">" + mainBondedTo[i].type + "!").c_str();
		}
		tmp = tmp + "*";

		return tmp;
	}
};
struct SpeciesExtended
{
	Coords mainAtom;
	std::vector <SpeciesExtended*> relations;
	char extraInfo_;

	SpeciesExtended(Coords crds)
	{
		mainAtom = crds;
		extraInfo_ = 'f';
	}

	SpeciesExtended&operator=(SpeciesExtended &specExt)
	{
		mainAtom = specExt.mainAtom;
		relations = specExt.relations;
		extraInfo_ = 'f';
		return *this;
	}
};



//Special form of the poscar class specifically for npa.  Inherits all other member functions / variables, and adds a few more
class npaPoscar : public Poscar
{
	public:
		//Variables
		std::vector <std::string> potSpeciesStr; ///to write to the top of the header file.  Values based on max nBonds and nAtom Types
		std::vector <Species> potSpeciesVect; ///vector of information provided in the string, except accessable
		std::vector <std::string> thisSpeciesStr; ///all existing species for this-> specific file, in string format
		std::vector <Species> thisSpeciesVect; ///all existing species for this-> file, as formatted as the Species struct
		std::vector <std::string> potDefectsStr;
		std::vector <SpeciesExtended*> thisSpeciesVectPtr;
		std::vector <std::vector <SpeciesExtended> > thisClustersVect; ///vector of existing clusters

		//Constructors, operators, etc
		npaPoscar();
		npaPoscar(std::string);
		npaPoscar(std::string, std::string);

		npaPoscar& operator= (const npaPoscar &npaposcar);

		//Function Declerations
		void fetchPotentialSpeciesString(int); ///generates the potential different relations between atoms, limited by the (int) max number of bonds on any atom.  Based off of the current input file being read in
		void fetchPotentialSpeciesString(int, std::string); ///same as above, except generates potential species based off of an input file
		void fetchPotentialSpeciesVector(); ///generates vector of species matching the potential species string
		void fetchPotentialSpeciesVector(std::string);
		void fetchThisSpeciesVector(std::string); ///generates species vector (for defects and non-defects on input file.  Based on existing bonds
		void fetchThisSpeciesVectorExtended(std::string);
		void fetchThisClustersVect(std::string); ///generates a vector of all clusters existing in model
		void writePotSpeciesStr(std::string, std::string); ///writes all potential species based on the potSpeciesStr
		void markPotSpeciesVect(std::string); ///given a user input file, marks the possible species as being defective or not (based on nBonds for a specific central atom)
		void writeThisFileBondDensity(std::string, std::string); ///writes specific file density in line with the potSpeciesStr - i.e., there will probably be a lot of 0s
		void writeThisFileBondDensity(std::string, std::string, std::string);
		void writeDefectInfo(std::string);
		void classifyPeriodicBonds(std::string);
		void clean(std::vector <std::string> &);
		std::vector <std::string> getRidOfUglyChars(std::vector <std::string>);
		bool periodicClusters(std::string); ///returns true if there is at least one instance of a periodic cluster (see definition of periodic cluster @ function definition)
		bool periodicClusters(std::vector <Species>); ///same as above, but checks for a specific group of species instead of the entire cell

		//stuff that works and shouldent be changed until you know what youre doing for sure
			void fetchThisSpeciesInfo(); ///generates the thisSpeciesStr and ThisSpeciesVect vectors 
			void fetchThisSpeciesInfo(std::string infileName); ///same as above, but does it with respect to an bondDataInfo's atoms instead of input file's atoms
};

//Constructor / Operator declerations
npaPoscar::npaPoscar() : Poscar()
{
	Poscar();
	potSpeciesStr;
	potSpeciesVect;
	thisSpeciesStr;
	thisSpeciesVect;
}

/*WHYYY is this not working???
npaPoscar::npaPoscar(const std::string path__) : Poscar(path__)
{
	Poscar(path__);
}*/

npaPoscar::npaPoscar(std::string inst, std::string path) : Poscar(inst, path)
{
	Poscar(inst, path);
	
	//Reads only the head of the file.  Does NOT return an iterator to the beginning of atomic coordinates; user will have to do this before reading those in
	if (inst == "ReadHead" || inst == "readHead")
	{
		//fetchPotentialSpeciesString(MAXBONDS);
		//potSpeciesVect;
	}

	//Reads entire file
	if (inst == "ReadAll" || inst == "readAll")
	{
		//fetchPotentialSpeciesString(MAXBONDS);
		//potSpeciesVect;
	}

}

npaPoscar& npaPoscar::operator=(const npaPoscar & npaposcar)
{
	Poscar::operator=(npaposcar);

	//Extra stuff here
	potSpeciesStr = npaposcar.potSpeciesStr;
	potSpeciesVect = npaposcar.potSpeciesVect;
	thisSpeciesStr = npaposcar.thisSpeciesStr;
	thisSpeciesVect = npaposcar.thisSpeciesVect;
	return *this;
}

//-------------------------------------------------------------------------------------------------

struct npaAtom
{
	std::string type;
	int count;
	std::string typeAndCount;

	void initTypeandCount()
	{
		std::stringstream newNameNum;
		newNameNum << count;
		typeAndCount = (newNameNum.str() + ">" + type + "!").c_str();
	}

	npaAtom()
	{
		type = "undf";
		count = 0;
	}

	npaAtom(std::string str, int cnt)
	{
		type = str;
		count = cnt;
		initTypeandCount();
	}
};

struct simpleSpecies ///for use with string initialization only
{
	std::string mainAtomType;
	std::vector <npaAtom> potentialBonds;

	simpleSpecies()
	{
		mainAtomType = "undf";
	}

	simpleSpecies(std::string str)
	{
		mainAtomType = str;
	}
};

//https://stackoverflow.com/questions/30290535/creating-n-nested-for-loops/30290814#30290814, my savior
void doThingWithNumber(const int* digits, int numDigits, int mod, std::vector <npaAtom> &vect, std::vector <std::string> atomTypes)
{
	int sum = 0;
	for (int i = numDigits - 1; i >= 0; i--)
		sum+= digits[i];

	if (sum < mod && sum > 0)
	{
		std::vector <int> digitsVect;
		for (int i = numDigits - 1, j = 0; i >= 0; i--, j++)
		{
			digitsVect.push_back(digits[i]);
		}
		for (int i = 0; i < digitsVect.size(); i++)
			vect.push_back(*new npaAtom(atomTypes[i], digitsVect[i]));
	}
}

void loopOverAllNumbers(int numDigits, int mod, std::vector <npaAtom> &vect, std::vector <std::string> atomTypes)
{
	int* digits = new int[numDigits];
	int i;
	for (int i = 0; i< numDigits; i++)
		digits[i] = 0;

	int maxDigit = 0;

	while (maxDigit < numDigits)
	{
		doThingWithNumber(digits, numDigits, mod, vect, atomTypes);
		for (i = 0; i < numDigits; i++) 
		{
			digits[i]++;
			if (digits[i] < mod)
				break;
			digits[i] = 0;
		}

		if (i > maxDigit)
			maxDigit = i;
	}
}

//initializes the vector holding all of the possible species for the poscar file that was read in, in string format
void npaPoscar::fetchPotentialSpeciesString(int maxBonds) ///assumes minBonds is 1
{	
	//If the string isnt empty, no need to run these calculations again (theoretically)
	if (!potSpeciesStr.empty())
		return;

	std::vector <std::vector <simpleSpecies> > mainAtoms;
	std::vector <simpleSpecies> potentialSpecies;
	std::vector <std::string> tmp = atomTypes;

	//sort atom types to avoid headaches
	std::sort(tmp.begin(), tmp.end(), npaStringSort);

	//Initialize species vector with names, not potential atoms for it to be bonded to
	for (int i = 0; i < atomTypes.size(); i++)
		potentialSpecies.push_back(*new simpleSpecies(atomTypes[i]));

	//For each element in the species vector, make a species list
	for (int i = 0; i < potentialSpecies.size(); i++)
	{
		loopOverAllNumbers(atomTypes.size(), maxBonds + 1, potentialSpecies[i].potentialBonds, tmp);
	}

	std::vector <std::string> speciesString;
	//make a vector of strings to print out
	for (int i = 0; i < atomTypes.size(); i++)
	{
		for (int j = 0; j < potentialSpecies[i].potentialBonds.size(); j = j + atomTypes.size()) ///this is sketchy; make sure to check for missed values
		{
			std::string tmp;
			for (int k = j; k < atomTypes.size() + j; k++)
			{
				tmp += potentialSpecies[i].potentialBonds[k].typeAndCount;
			}
			speciesString.push_back(*new std::string = (atomTypes[i] + ":" + tmp + "*").c_str());

		}
	}

	//Delete all instances of H being bonded to more than one atom
	int sum = 0;
	for (int i = 0; i < speciesString.size(); i++)
	{
		if (speciesString[i][0] == 'H')
			for (int j = 0; j < speciesString[i].size(); j++)
				if (speciesString[i][j] == ':' || speciesString[i][j] == '!')
					if (speciesString[i][j+1] != '0' && speciesString[i][j + 1] != '*')
						sum += speciesString[i][j+1] - '0'; ///if you get errors, you'll need a more sophisticated way of determining the sum.  This should suffice for NPA though
		if (sum > 1)										///also, I forsee a possible error with the conversion from char to int.  May be worth casting the char as an int before math
			speciesString[i] = "del";
		sum = 0;
	}

	std::vector <std::string> tmp2;
	for (int i = 0; i < speciesString.size(); i++)
		if (speciesString[i] != "del")
			tmp2.push_back(speciesString[i]);

	speciesString = tmp2;

	std::sort(speciesString.begin(), speciesString.end(), npaStringSort);
	potSpeciesStr = speciesString;
}

void npaPoscar::fetchPotentialSpeciesString(int maxBonds, std::string inFileName)
{
	//If the string isnt empty, no need to run these calculations again (theoretically)
	if (!potSpeciesStr.empty())
		return;

	//I just copied and pasted this function from the one above, so to make that work, i'm reading in the atoms from input file and pretending
	//that theyre the atoms in an npaPOSCAR file head
	std::ifstream infile(inFileName.c_str()); ///infileName should be the bond data info, by the way
		if (infile.fail()){
			std::cout << "fetchPotentialSpeciesString: could not open input file for reading\n";
			return;}

	//Read in the atom types from infilename
	std::vector <std::string> bondInfoNames;

        //I shouldent have to do this
        int nIt = 0;

		std::ifstream infile_(inFileName.c_str());
		while(std::getline(infile_, *new std::string))
			nIt++;
		infile_.close();
        
        for (int i = 0; i < nIt; i++)
        {
        	std::string str1, str2, garbage;
        
                infile >> str1 >> str2 >> garbage;
                bondInfoNames.push_back(str1);
                bondInfoNames.push_back(str2);
         }
        
	std::vector <std::string> tempo;
	for (int i = 0; i < bondInfoNames.size(); i++)
		for (int j = 0; j < bondInfoNames.size(); j++)
			if (i != j)
				if (bondInfoNames[i] == bondInfoNames[j])
					bondInfoNames[j] = "del";
	for (int i = 0; i < bondInfoNames.size(); i++)
		if (bondInfoNames[i] != "del")
			tempo.push_back(bondInfoNames[i]);

	bondInfoNames = tempo;

	//make a npaPoscar instance, and call set the atomTypes to the atoms read in from the infile
	npaPoscar tempNPA;
	tempNPA.atomTypes = bondInfoNames;

	//Now, the rest is exactly the same as the above function
	std::vector <std::vector <simpleSpecies> > mainAtoms;
	std::vector <simpleSpecies> potentialSpecies;
	std::vector <std::string> tmp = tempNPA.atomTypes;

	//sort atom types to avoid headaches
	std::sort(tmp.begin(), tmp.end(), npaStringSort);

	//Initialize species vector with names, not potential atoms for it to be bonded to
	for (int i = 0; i < tempNPA.atomTypes.size(); i++)
		potentialSpecies.push_back(*new simpleSpecies((tempNPA.atomTypes[i]).c_str()));

	//For each element in the species vector, make a species list
	for (int i = 0; i < potentialSpecies.size(); i++)
	{
		loopOverAllNumbers(tempNPA.atomTypes.size(), maxBonds + 1, potentialSpecies[i].potentialBonds, tmp);
	}

	std::vector <std::string> speciesString;
	//make a vector of strings to print out
	for (int i = 0; i < tempNPA.atomTypes.size(); i++)
	{
		for (int j = 0; j < potentialSpecies[i].potentialBonds.size(); j = j + tempNPA.atomTypes.size()) ///this is sketchy; make sure to check for missed values
		{
			std::string tmp;
			for (int k = j; k < tempNPA.atomTypes.size() + j; k++)
			{
				tmp += potentialSpecies[i].potentialBonds[k].typeAndCount;
			}
			speciesString.push_back(*new std::string = (tempNPA.atomTypes[i] + ":" + tmp + "*").c_str());

		}
	}

	//Delete all instances of bonds that are defective (current:  any C, Si != 4 & H != 1)
	//TODO: fix this - it is very lazy and inefficient and not generalized at all
	//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	int sum = 0;
	for (int i = 0; i < speciesString.size(); i++)
	{
		if (speciesString[i][0] == 'H')
			for (int j = 0; j < speciesString[i].size(); j++)
				if (speciesString[i][j] == ':' || speciesString[i][j] == '!')
					if (speciesString[i][j + 1] != '0' && speciesString[i][j + 1] != '*')
						sum += speciesString[i][j + 1] - '0'; ///if you get errors, you'll need a more sophisticated way of determining the sum.  This should suffice for NPA though
		if (sum > 1)										///also, I forsee a possible error with the conversion from char to int.  May be worth casting the char as an int before math
		{	
			potDefectsStr.push_back(speciesString[i]);
			speciesString[i] = "del";
		}
		sum = 0;
	}

	std::vector <std::string> tmp2;
	for (int i = 0; i < speciesString.size(); i++)
		if (speciesString[i] != "del")
			tmp2.push_back(speciesString[i]);

	speciesString = tmp2;
	
	sum = 0;
	for (int i = 0; i < speciesString.size(); i++)
	{
		if (speciesString[i][0] == 'C')
			for (int j = 0; j < speciesString[i].size(); j++)
				if (speciesString[i][j] == ':' || speciesString[i][j] == '!')
					if (speciesString[i][j + 1] != '0' && speciesString[i][j + 1] != '*')
						sum += speciesString[i][j + 1] - '0'; ///if you get errors, you'll need a more sophisticated way of determining the sum.  This should suffice for NPA though
		if (sum != 4 && speciesString[i][0] == 'C')										///also, I forsee a possible error with the conversion from char to int.  May be worth casting the char as an int before math
		{
			potDefectsStr.push_back(speciesString[i]);
			speciesString[i] = "del";
		}
		sum = 0;
	}

	std::vector <std::string> tmp3;
	for (int i = 0; i < speciesString.size(); i++)
		if (speciesString[i] != "del")
			tmp3.push_back(speciesString[i]);

	speciesString = tmp3;

	sum = 0;
	for (int i = 0; i < speciesString.size(); i++)
	{
		if (speciesString[i][0] == 'S')
			for (int j = 0; j < speciesString[i].size(); j++)
				if (speciesString[i][j] == ':' || speciesString[i][j] == '!')
					if (speciesString[i][j + 1] != '0' && speciesString[i][j + 1] != '*')
						sum += speciesString[i][j + 1] - '0'; ///if you get errors, you'll need a more sophisticated way of determining the sum.  This should suffice for NPA though
		if (sum != 4 && speciesString[i][0] == 'S')										///also, I forsee a possible error with the conversion from char to int.  May be worth casting the char as an int before math
		{
			potDefectsStr.push_back(speciesString[i]);
			speciesString[i] = "del";
		}
		sum = 0;
	}

	std::vector <std::string> tmp4;
	for (int i = 0; i < speciesString.size(); i++)
		if (speciesString[i] != "del")
			tmp4.push_back(speciesString[i]);

	speciesString = tmp4;

	sum = 0;
	for (int i = 0; i < speciesString.size(); i++)
		if (speciesString[i][0] == 'H')
			speciesString[i] = "del";

	std::vector <std::string> tmp5;
	for (int i = 0; i < speciesString.size(); i++)
		if (speciesString[i] != "del")
			tmp5.push_back(speciesString[i]);

	speciesString = tmp5;

	sum = 0;
	for (int i = 0; i < speciesString.size(); i++)
	{
		if (speciesString[i][0] == 'C' || speciesString[i][0] == 'S')
			for (int j = 0; j < speciesString[i].size(); j++)
				if (speciesString[i][j] == 'H' && speciesString[i][j-2] != '0')
						sum += speciesString[i][j - 2] - '0'; ///if you get errors, you'll need a more sophisticated way of determining the sum.  This should suffice for NPA though
		if (sum >= 4 && (speciesString[i][0] == 'C' || speciesString[i][0] == 'S'))										///also, I forsee a possible error with the conversion from char to int.  May be worth casting the char as an int before math
		{
			potDefectsStr.push_back(speciesString[i]);
			speciesString[i] = "del";
		}
		sum = 0;
	}

	std::vector <std::string> tmp6;
	for (int i = 0; i < speciesString.size(); i++)
		if (speciesString[i] != "del")
			tmp6.push_back(speciesString[i]);

	speciesString = tmp6;
	
	//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	std::sort(speciesString.begin(), speciesString.end(), npaStringSort);
	std::sort(potDefectsStr.begin(), potDefectsStr.end(), npaStringSort);
	potSpeciesStr = speciesString;

	infile.close();
}

//initializes the vector of accessable species info
void npaPoscar::fetchPotentialSpeciesVector()
{
	if (!potSpeciesVect.empty())
		return;

	fetchPotentialSpeciesString(MAXBONDS);
	for (int i = 0; i < potSpeciesStr.size(); i++)
		potSpeciesVect.push_back(*new Species ((potSpeciesStr[i]).c_str()));
}

//same as above, but calls the string version
void npaPoscar::fetchPotentialSpeciesVector(std::string inFileName)
{
	if (!potSpeciesVect.empty())
		return;

	fetchPotentialSpeciesString(MAXBONDS, inFileName.c_str());
	for (int i = 0; i < potSpeciesStr.size(); i++)
		potSpeciesVect.push_back(*new Species((potSpeciesStr[i]).c_str()));
}

void npaPoscar::fetchThisSpeciesVector(std::string infile)
{
	if (!thisSpeciesVect.empty())
		return;

	fetchAtomBonds(infile.c_str());

	//fill main atoms
	for (int i = 0; i < atomCoords.size(); i++)
	{
		thisSpeciesVect.push_back(*new Species(atomCoords[i]));
		thisSpeciesVect[i].mainAtom.extraInfo = "central";
	}

	//fill atom bonds, with the main atom ALWAYS the first atom listed in vector of atom bonds
	for (int i = 0; i < thisSpeciesVect.size(); i++)
		for (int j = 0; j < atomBonds.size(); j++)
			if (thisSpeciesVect[i].mainAtom == atomBonds[j].pairedAtoms[0])
				thisSpeciesVect[i].atomBonds_.push_back(atomBonds[j]);
			else
				if (thisSpeciesVect[i].mainAtom == atomBonds[j].pairedAtoms[1])
				{
					atomPair tmp = atomBonds[j];
					Coords tmpCrds = tmp.pairedAtoms[1];
					tmp.pairedAtoms[1] = tmp.pairedAtoms[0];
					tmp.pairedAtoms[0] = tmpCrds;
					thisSpeciesVect[i].atomBonds_.push_back(tmp);
				}

}

void clusterCounter(SpeciesExtended &specExt, std::vector <SpeciesExtended> &vect)
{
	
}
/*
void npaPoscar::fetchThisSpeciesVectorExtended(std::string infile)
{
	if (!thisSpeciesVectPtr.empty())
		return;

	fetchThisSpeciesVector(infile.c_str());
	std::vector <SpeciesExtended*> extSpeciesVect;

	//Initialize main atoms and their relations to other atoms
	for (int i = 0; i < thisSpeciesVect.size(); i++)
	{
		extSpeciesVect.push_back(new SpeciesExtended(*new Coords(thisSpeciesVect[i].mainAtom)));
		for (int j = 0; j < thisSpeciesVect[i].atomBonds_.size(); j++)
			extSpeciesVect[i]->relations.push_back(new SpeciesExtended(*new Coords (thisSpeciesVect[i].atomBonds_[j].pairedAtoms[1])));
	}

	std::cout << "";

	std::cout << "";

	//Link together all atoms through pointers
	for (int i = 0; i < extSpeciesVect.size(); i++)
		for (int j = 0; j < extSpeciesVect.size(); j++)
			for (int k = 0; k < extSpeciesVect[j]->relations.size(); k++)
			{
				if (extSpeciesVect[i]->mainAtom == extSpeciesVect[j]->relations[k]->mainAtom)
				{
					*extSpeciesVect[j]->relations[k] = *extSpeciesVect[i];
				}
			}

	std::cout << "";

	std::cout << "";

	thisSpeciesVectPtr = extSpeciesVect;
	extSpeciesVect[0]->extraInfo_ = 'y';
	thisSpeciesVectPtr[0]->extraInfo_ = 'y';

	std::cout << "";

	std::cout << "";
}

void DFS_(SpeciesExtended &spec, std::vector <SpeciesExtended> &vect)
{
	spec.extraInfo_ = 'i'; ///i for 'intermediate'

	for (int i = 0; i < spec.relations.size(); i++)
		if (spec.relations[i]->extraInfo_ == 'u') ///u for 'uncounted'
		{
			vect.push_back(*spec.relations[i]);
			DFS_(*spec.relations[i], vect);
		}

	spec.extraInfo_ = 'd'; ///d for 'done'
	return;
}

void npaPoscar::fetchThisClustersVect(std::string infile)
{
	if (!thisClustersVect.empty())
		return;

	fetchThisSpeciesVectorExtended(infile.c_str());

	//Set all atoms to uncounted
	for (int i = 0; i < thisSpeciesVectPtr.size(); i++)
		thisSpeciesVectPtr[i]->extraInfo_ = 'u';
	
	//Determine Clusters, fill vector
	for (int i = 0; i < thisSpeciesVectPtr.size(); i++)
	{
		std::vector <SpeciesExtended> tmpCluster;
		if (thisSpeciesVectPtr[i]->extraInfo_ == 'u')
		{
			tmpCluster.push_back(*thisSpeciesVectPtr[i]);
			DFS_(*thisSpeciesVectPtr[i], tmpCluster);
		}
		thisClustersVect.push_back(tmpCluster);
	}

	std::cout << "";

	std::cout << "";

}
*/
void npaPoscar::classifyPeriodicBonds(std::string infile)
{
	fetchAtomBonds(infile.c_str());

	for (int i = 0; i < atomBonds.size(); i++)
		for (int n = 0; n < 2; n++)
			if (atomBonds[i].pairedAtoms[n].extraInfo != "original")
				atomBonds[i].extraInfo = "periodic";
			else
				atomBonds[i].extraInfo = "inCell";
}


struct coordsAndPartners
{
	Coords thisAtom;
	std::vector <Coords*> thisAtomPartners;
	int algProgress = 1; ///number to indicate if an atom needs to be considered or not

	coordsAndPartners()
	{
		Coords thisAtom;
	}

	coordsAndPartners(Coords crds)
	{
		thisAtom = crds;
	}
};


void npaPoscar::fetchThisSpeciesInfo()
{
	if (!thisSpeciesVect.empty())
		return;

	fetchAtomBonds();
	std::vector <Species> allAtoms;

	for (int i = 0; i < atomCoords.size(); i++)
		allAtoms.push_back(*new Species(atomCoords[i]));

	//Filling the bonds vector for each atom species in allAtoms vector
	for (int i = 0; i < allAtoms.size(); i++)
		for (int j = 0; j < atomBonds.size(); j++)
				if (allAtoms[i].mainAtom == atomBonds[j].pairedAtoms[0] || allAtoms[i].mainAtom == atomBonds[j].pairedAtoms[1])
					allAtoms[i].atomBonds_.push_back(atomBonds[j]);
	
	thisSpeciesVect = allAtoms; ///TODO:  'allAtoms' isnt even necessary, could just initialize thisSpeciesVect right away without making a new vector
								///when you have time, get rid of it to save memory

	if (!thisSpeciesStr.empty())
		return;

	std::vector <std::string> toWrite;
	//Filling a string vector with the string versions of the found species
	for (int i = 0; i < allAtoms.size(); i++)
		thisSpeciesStr.push_back(allAtoms[i].writeString(atomTypes));
}

void npaPoscar::fetchThisSpeciesInfo(std::string inFileName)
{
	if (!thisSpeciesVect.empty())
		return;

	//I just copied and pasted this function from the one above, so to make that work, i'm reading in the atoms from input file and pretending
	//that theyre the atoms in an npaPOSCAR file head
	std::ifstream infile(inFileName.c_str()); ///infileName should be the bond data info, by the way
	if (infile.fail()){
		std::cout << "fetchThisSpeciesInfo: could not open input file for reading\n";
		return;}

	//Read in the atom types from infilename
	std::vector <std::string> bondInfoNames;

 //I shouldent have to do this but g++ is yet again making me do some crazy workaround just to properly read in files
       int nIt = 0;
 
       std::ifstream infile_(inFileName.c_str());
       while(std::getline(infile_, *new std::string))
       		nIt++;
       infile_.close();

      for (int i = 0; i < nIt; i++)
	  {
		std::string str1, str2, garbage;

		infile >> str1 >> str2 >> garbage;
		bondInfoNames.push_back(str1);
        bondInfoNames.push_back(str2);
      }
 
	std::vector <std::string> tempo;
	for (int i = 0; i < bondInfoNames.size(); i++)
		for (int j = 0; j < bondInfoNames.size(); j++)
			if (i != j)
				if (bondInfoNames[i] == bondInfoNames[j])
					bondInfoNames[j] = "del";
	for (int i = 0; i < bondInfoNames.size(); i++)
		if (bondInfoNames[i] != "del")
			tempo.push_back(bondInfoNames[i]);
	std::sort(tempo.begin(), tempo.end(), npaStringSort);

	bondInfoNames = tempo;
	//Create new instance of npaCoords, and change the atomTypes
	npaPoscar npaCOORDS = *this;
	npaCOORDS.atomTypes = bondInfoNames;

	npaCOORDS.fetchAtomBonds(inFileName.c_str());

	std::vector <Species> allAtoms;

	
	for (int i = 0; i < npaCOORDS.atomCoords.size(); i++)
		allAtoms.push_back(*new Species(npaCOORDS.atomCoords[i]));

	//Filling the bonds vector for each atom species in allAtoms vector
	for (int i = 0; i < allAtoms.size(); i++)
		for (int j = 0; j < npaCOORDS.atomBonds.size(); j++)
			if (allAtoms[i].mainAtom == npaCOORDS.atomBonds[j].pairedAtoms[0] || allAtoms[i].mainAtom == npaCOORDS.atomBonds[j].pairedAtoms[1])
				allAtoms[i].atomBonds_.push_back(npaCOORDS.atomBonds[j]);

	thisSpeciesVect = allAtoms; ///TODO:  'allAtoms' isnt even necessary, could just initialize thisSpeciesVect right away without making a new vector
								///when you have time, get rid of it

	if (!thisSpeciesStr.empty())
		return;

	std::vector <std::string> toWrite;

	//Filling a string vector with the string versions of the found species
	for (int i = 0; i < allAtoms.size(); i++)
		thisSpeciesStr.push_back(allAtoms[i].writeString(npaCOORDS.atomTypes));
}

//marks defective and nondefective atoms in potSpeciesVect
void npaPoscar::markPotSpeciesVect(std::string infilePath)
{
	std::ifstream infile(infilePath.c_str());
	if (infile.fail())
	{
		std::cout << "markPotSpeciesVect: failed to open " << infilePath << "\n";
		return;
	}
}

//writes potential species to (normally) the top of the csv file
void npaPoscar::writePotSpeciesStr(std::string outfilePath, std::string bondDataInfoLoc)
{
	fetchPotentialSpeciesString(MAXBONDS, bondDataInfoLoc.c_str());

	std::ofstream outfile((outfilePath + ".csv").c_str());
	std::ofstream outfileDefect((outfilePath + "defects" + ".csv").c_str());

	std::sort(potSpeciesStr.begin(), potSpeciesStr.end(), npaStringSort);
	std::sort(potDefectsStr.begin(), potDefectsStr.end(), npaStringSort);

	std::vector <std::string> cleaned = getRidOfUglyChars(potSpeciesStr);
	std::vector <std::string> cleanedDefects = getRidOfUglyChars(potDefectsStr);

	//write to file
	for (int i = 0; i < cleaned.size(); i++)
		outfile << "," << cleaned[i];

	//write to file
	for (int i = 0; i < cleanedDefects.size(); i++)
		outfileDefect << "," << cleanedDefects[i];

	outfile.close();
	outfileDefect.close();
}

//writes bond density to the csv file
void npaPoscar::writeThisFileBondDensity(std::string instructions, std::string outfilePath)
{
	if (instructions[0] != 'w')
	{
		std::cout << "writeThisFileBondDensity: the first argument to this function should be the instructions:  writeAll or writeDefects\n";
		return;
	}

	//If user wants to write all possible types
	if (instructions == "writeAll" || instructions == "writeall")
	{
		fetchThisSpeciesInfo();

		std::ofstream outfile(outfilePath.c_str());
		std::ofstream outfileDefect((outfilePath + "defects").c_str());

		std::sort(potSpeciesStr.begin(), potSpeciesStr.end(), npaStringSort);
		std::sort(potDefectsStr.begin(), potDefectsStr.end(), npaStringSort);

		//Will be using the previously defined atomTypesandCounts struct to keep track of the counts easily
		std::vector <atomTypesAndCounts> speciesTypesAndCounts;
		for (int i = 0; i < potSpeciesStr.size(); i++)
			speciesTypesAndCounts.push_back(*new atomTypesAndCounts((potSpeciesStr[i]).c_str()));

		//Initialize counts
		for (int i = 0; i < speciesTypesAndCounts.size(); i++)
			for (int j = 0; j < thisSpeciesStr.size();j++)
				if (speciesTypesAndCounts[i].type == thisSpeciesStr[j])
					speciesTypesAndCounts[i].count++;

		//write to file
		outfile << "\n";
		for (int i = 0; i < speciesTypesAndCounts.size(); i++)
			outfile << "," << speciesTypesAndCounts[i].count;

		outfile.close();
	}
}

void npaPoscar::writeThisFileBondDensity(std::string instructions, std::string outfilePath, std::string infilePath)
{
	if (instructions[0] != 'w')
	{
		std::cout << "writeThisFileBondDensity: the first argument to this function should be the instructions:  writeAll or writeDefects\n";
		return;
	}

	//If user wants to write all possible types
	if (instructions == "writeAll" || instructions == "writeall")
	{
		fetchThisSpeciesInfo(infilePath.c_str()); ///open bondDataInfo
		std::ofstream outfile((outfilePath + ".csv").c_str());
		std::ofstream outfileDefect((outfilePath + "defects" + ".csv").c_str());

		std::sort(potSpeciesStr.begin(), potSpeciesStr.end(), npaStringSort);
		std::sort(potDefectsStr.begin(), potDefectsStr.end(), npaStringSort);

		//Clean
		clean(potSpeciesStr);
		clean(potDefectsStr);

		//Will be using the previously defined atomTypesandCounts struct to keep track of the counts easily
		std::vector <atomTypesAndCounts> speciesTypesAndCounts;
		std::vector <atomTypesAndCounts> speciesTypesAndCountsDefects;
		for (int i = 0; i < potSpeciesStr.size(); i++)
			speciesTypesAndCounts.push_back(*new atomTypesAndCounts((potSpeciesStr[i]).c_str()));
		for (int i = 0; i < potDefectsStr.size(); i++)
			speciesTypesAndCountsDefects.push_back(*new atomTypesAndCounts((potDefectsStr[i]).c_str()));

		//Initialize counts
		for (int i = 0; i < speciesTypesAndCounts.size(); i++)
			for (int j = 0; j < thisSpeciesStr.size(); j++)
				if (speciesTypesAndCounts[i].type == thisSpeciesStr[j])
					speciesTypesAndCounts[i].count++;

		//Initialize counts
		for (int i = 0; i < speciesTypesAndCountsDefects.size(); i++)
			for (int j = 0; j < thisSpeciesStr.size(); j++)
				if (speciesTypesAndCountsDefects[i].type == thisSpeciesStr[j])
					speciesTypesAndCountsDefects[i].count++;

		//write to file
		outfile << "\n";
		for (int i = 0; i < speciesTypesAndCounts.size(); i++)
			outfile << "," << speciesTypesAndCounts[i].count;

		//write to file
		outfileDefect << "\n";
		for (int i = 0; i < speciesTypesAndCountsDefects.size(); i++)
			outfileDefect << "," << speciesTypesAndCountsDefects[i].count;

		outfile.close();
		outfileDefect.close();
	}
}

void npaPoscar::clean(std::vector <std::string> & vectToClean)
{
	for (int i = 0; i < vectToClean.size(); i++)
		if (vectToClean[i][0] == ':')
			vectToClean[i] == "rem";

	std::vector <std::string> tmp;
	for (int i = 0; i < vectToClean.size(); i++)
		if (vectToClean[i] != "rem")
			tmp.push_back(vectToClean[i]);

	vectToClean = tmp;
}

std::vector <std::string> npaPoscar::getRidOfUglyChars(std::vector <std::string> vectToClean)
{
	std::vector <std::string> toReturn;

	for (int i = 0; i < vectToClean.size(); i++)
	{
		std::vector <char> charVect;
		for (int j = 0; j < vectToClean[i].size(); j++)
		{
			
			if (isalnum(vectToClean[i][j]) || vectToClean[i][j] == ':')
				charVect.push_back(vectToClean[i][j]);
		}
		toReturn.push_back((*new std::string(charVect.begin(), charVect.end())));
	}


	if (toReturn.size() != vectToClean.size())
		std::cout << "getRidOfUglyChars:  something awful has happened\n";

	return toReturn;
}

//Checks all legitimate bonds for an instance of two bonded atoms being the same type, but one of them being an "original" atom and the other being
//anything but an original atom - Original atoms are the ones that are in the supercell, as defined in the 'extend supercell' function
//Returns true if at least one A-A bond is found where either A1 or A2 is original, and the other is not
///IMPORTANT NOTE:  this relies on the 'extraInfo' section of the Coords struct to be unaltered since the last calling of 'expantSupercell()'.  
///make sure to use this check before altering the 'extraInfo' section of any atom
bool npaPoscar::periodicClusters(std::string bondInfoLoc)
{
	bool periodicCluster = false;
	fetchAtomBonds(bondInfoLoc.c_str());

	for (int i = 0; i < atomBonds.size(); i++)
		if (atomBonds[i].pairedAtoms[0].atomType == atomBonds[i].pairedAtoms[1].atomType) ///if the atom types are the same (necessary for something to be considered a cluster)
			if (atomBonds[i].pairedAtoms[0].extraInfo != "original" || atomBonds[i].pairedAtoms[1].extraInfo != "original") ///then if any of the two bonded atoms are non-original...
				if (atomBonds[i].pairedAtoms[0].extraInfo != "original" || atomBonds[i].pairedAtoms[1].extraInfo != "original") ///and any of the other two bonded atoms are original
					if (atomBonds[i].pairedAtoms[0].extraInfo != atomBonds[i].pairedAtoms[1].extraInfo) ///and the two types arent equal (they shouldent be at this point, but just adding this line in for saftey against false positives
						periodicCluster = true; ///then there must be a bond across at least one boundry

	return periodicCluster;
}

//Same info (and warnings) that were given for te previous periodicCluster check apply to this one
bool npaPoscar::periodicClusters(std::vector <Species> specVect)
{
	bool periodicCluster = false;

	for (int i = 0; i < specVect.size(); i++)
		for (int j = 0; j < specVect[i].atomBonds_.size(); j++)
			if (specVect[i].atomBonds_[j].pairedAtoms[0].atomType == specVect[i].atomBonds_[j].pairedAtoms[1].atomType) ///if the atom types are the same (necessary for something to be considered a cluster)
				if (specVect[i].atomBonds_[j].pairedAtoms[0].extraInfo != "original" || specVect[i].atomBonds_[j].pairedAtoms[1].extraInfo != "original") ///then if any of the two bonded atoms are non-original...
					if (specVect[i].atomBonds_[j].pairedAtoms[0].extraInfo != "original" || specVect[i].atomBonds_[j].pairedAtoms[1].extraInfo != "original") ///and any of the other two bonded atoms are original
						if (specVect[i].atomBonds_[j].pairedAtoms[0].extraInfo != specVect[i].atomBonds_[j].pairedAtoms[1].extraInfo) ///and the two types arent equal (they shouldent be at this point, but just adding this line in for saftey against false positives
							periodicCluster = true; ///then there must be a bond across at least one boundry

	return periodicCluster;
}
