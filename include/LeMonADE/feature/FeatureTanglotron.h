/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by 
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers (see AUTHORS)
    ooo                        | 
----------------------------------------------------------------------------------

This file is part of LeMonADE.

LeMonADE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LeMonADE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LeMonADE.  If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------------*/

#ifndef LEMONADE_FEATURE_TANGLOTRON_H
#define LEMONADE_FEATURE_TANGLOTRON_H

/**
 * @file
 * @date 2021/02/03
 * @author Martin Wengenmayr/Cornelia Struebig
 * @brief Definition and Implementation of FeatureTanglotron
**/
#include <iostream>
#include <cmath>

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/feature/FeatureBoltzmann.h>
#include <LeMonADE/feature/FeatureTanglotronReadWrite.h>

/**
* @class TanglotronAttributeTag
* 
* @brief Extension for tanglotrons - isTanglotron and tanglotronType
* 
* @details Extends tanglotrons by a bool as tag isTanglotron along with getter and
* setter. Initially this tag is set to false.
* Besides another enum tag tanglotronType is set, which holds information about
* stator (int 0) or rotorA,B,C,D (int 1-4). Also this goes along with getter and setter
* functions. This tag is initially set to 5.
**/
class TanglotronAttributeTag
{
public:
    //! enum variable holding the tanglotronType
    enum tanglotronType{
        stator=0, 
        rotorA=1, 
        rotorB=2,
        rotorC=3,
        rotorD=4
    };
   
	//! Standard constructor- initially the tag is set to 5 and isTanglotron is set to false
	TanglotronAttributeTag():isTanglotron(false),tag(5),tanglotronID(){}

	/**
	* @brief setting boolean tag isTanglotron
	*
	* @param isTanglotron_ bool which tells if monomer is part of a tanglotron
    * motor or not
	**/
	void setIsTanglotron(bool isTanglotron_){ isTanglotron=isTanglotron_; }
	
    //! getting boolean tag isTanglotron
    bool getIsTanglotron() const {return isTanglotron;}
    
    /**
	* @brief setting enum tag tanglotronType
	*
	* @param attribute role of tanglotron part as string (stator, rotorA,B,C,D) or
    * integer (0,1,2,3,4)
	**/
    void setTanglotronType(int attribute){ tag=attribute; }
    
	//! getting enum tag tanglotronType
	int32_t getTanglotronType() const {if(isTanglotron) return tag;}
	
	/**
	* @brief setting size_t tag tanglotronID
	*
	* @param id number of motor in vector TanglotronMotorUnits
	**/
    void setTanglotronID(size_t id){ tanglotronID=id; }
    
    //! getting size_t tag tanglotronID
	size_t getTanglotronID() const {if(isTanglotron) return tanglotronID;}
	
private:
    //! Private Variable holding the bool if monomer is part of tanglotron
    bool isTanglotron;
	//! Private variable holding the tanglotronType. Default is 5.
    int32_t tag;
    //! private variable holding the number of one tanglotron motor in the vector holding all motors
    size_t tanglotronID;
};

/**
* @class TanglotronMotor
*
* @brief set and get indices of stator and rotors of a tanglotron motor
* 
* @details One tanglotron motor is built out of one central stator monomer and
* four surrounding rotor monomers which are connected each with the stator
* and in pairs with each other (A and B as well as C and D).
*
**/
class TanglotronMotor: public Feature
{
public:
    
    //! Standard constructor- initially all entries are set to NULL
	TanglotronMotor():StatorIndex(0), RotorAIndex(0), RotorBIndex(0), RotorCIndex(0), RotorDIndex(0), LK(0), LK_progress(), angleSum(0.0), angle_progress(), angleSumAB(0.0), angleSumCD(0.0) {}

    //! setter for tanglotron stator/rotors
	void setTanglotronMotor(uint32_t idxStator, uint32_t idxA, uint32_t idxB, uint32_t idxC, uint32_t idxD);

    //! getter for StatorIndex
    uint32_t getIndexStator() const{return StatorIndex;}
    //! getter for RotorAIndex
    uint32_t getIndexRotorA() const{return RotorAIndex;}
    //! getter for RotorBIndex
    uint32_t getIndexRotorB() const{return RotorBIndex;}
    //! getter for RotorCIndex
	uint32_t getIndexRotorC() const{return RotorCIndex;}
    //! getter for RotorDIndex
    uint32_t getIndexRotorD() const{return RotorDIndex;}

    //! getter for LK
    int64_t getLK() const{return LK;}
    //! getter for LK_progress
    std::vector<int8_t> getLKprogress() const{return LK_progress;}

    //! getter for angleSum
    double getAngleSum() const{return angleSum;}
    //! getter for angle_progress
    std::vector<double> getAngleSumProgress() const{return angle_progress;}
    //! getter for angleSumAB
    double getAngleSumAB() const{return angleSumAB;}
    //! getter for angleSumCD
    double getAngleSumCD() const{return angleSumCD;}
    
    //! update counter of linking number
    void addLK(int8_t halfLKtoAddOrSubstract);
    //! update full tanglotron angle
    void addAngle(double angle);
    //! update angle stator/AB rotor unit
    void addAngleAB(double angle){angleSumAB += angle;}
    //! update angle stator/CD rotor unit
    void addAngleCD(double angle){angleSumCD += angle;}
    
private:
    
    //! Index of the stator monomer of the tanglotron motor
    uint32_t StatorIndex;
    //! Index of the first (A) rotor monomer of the tanglotron motor
    uint32_t RotorAIndex;
    //! Index of the second (B) rotor monomer of the tanglotron motor
    uint32_t RotorBIndex;
    //! Index of the third (C) rotor monomer of the tanglotron motor
    uint32_t RotorCIndex;
    //! Index of the fourth (D) rotor monomer of the tanglotron motor
    uint32_t RotorDIndex;
    
    //! Linking Number, which will be counted on the run
    int64_t LK;
    //! Vector, which holds all added or substracted half Linking Numbers stepwise
    std::vector<int8_t> LK_progress;
    
    //! Linking Number, which will be counted on the run
    double angleSum;
    //! Vector, which holds current angleSum stepwise
    std::vector<double> angle_progress;
    //! motor angle, which will be counted on the run
    double angleSumAB;
    //! motor angle, which will be counted on the run
    double angleSumCD;
};

/**
* @details function to set the indices of the tanglotron motor unit
*/
void TanglotronMotor::setTanglotronMotor(uint32_t idxStator, uint32_t idxA, uint32_t idxB, uint32_t idxC, uint32_t idxD)
{
	StatorIndex = idxStator;
	RotorAIndex = idxA;
	RotorBIndex = idxB;
    RotorCIndex = idxC;
    RotorDIndex = idxD;
}

/**
* @details function to add or substract half Linking Numbers
*/
void TanglotronMotor::addLK(int8_t halfLKtoAddOrSubstract)
{
    #ifdef DEBUG
        std::cout << "in addLK with " << halfLKtoAddOrSubstract << std::endl;
    #endif /*DEBUG*/
	LK += halfLKtoAddOrSubstract;
    LK_progress.push_back(halfLKtoAddOrSubstract);
}

/**
* @details function to add or substract angle from angleSum
*/
void TanglotronMotor::addAngle(double angle)
{
	angleSum += angle;
    angle_progress.push_back(angle);
}


/**
* @class FeatureTanglotron
*
* @brief apply potential for one directional rotation of tanglotron motor
* 
* @details If one of four rotor monomers of the tanglotron is moved by simulation
* Metropolis-criterion is calculated by torque and angle of the move.
*
* @tparam IngredientsType
**/
class FeatureTanglotron: public Feature
{
public:
    //! requested feature: FeatureBoltzmann to calculate move probability from calculated potential energy
    typedef LOKI_TYPELIST_1(FeatureBoltzmann) required_features_back;
    
    //! monomer_extensions, gives information if monomer is part of tanglotron motors and special which part
    typedef LOKI_TYPELIST_1(TanglotronAttributeTag) monomer_extensions;
    
    /**
    * constructor
    *
    * @param variableCalculatePotential to decide if tanglotron potential should
    * be calculated or not (default true)
    */
    FeatureTanglotron():variableCalculatePotential(true) {}
    
    //! default empty destructor
    virtual ~FeatureTanglotron(){}
    
    /**
    * synchronize
    *
    * @param ingredients A reference to the IngredientsType - mainly the system.
    */
    template<class IngredientsType> void synchronize(IngredientsType& ingredients)
    {
        #ifdef DEBUG
            std::cout << "Number of Tanglotrons in Vector: " << getTanglotronMotors().size() << std::endl;
        #endif /*DEBUG*/
    }
    
    //! check move with calculation of tanglotron potential energy
	template<class IngredientsType> 
	bool checkMove(const IngredientsType& ingredients,MoveLocalSc& move);
    
    //! default behavior for unknown moves: return true
    template<class IngredientsType> 
	bool checkMove(const IngredientsType& ingredients,MoveBase& move){ return true; }
	
	//! getter for variableCalculatePotential
	bool getCalculatePotential() const {return variableCalculatePotential;}
	
	//! setter for variableCalculatePotential
    void setCalculatePotential(bool yesorno){variableCalculatePotential=yesorno;}
    
    //! getter for Torque
    double getTorque() const{return Torque;}
    
    //! setter for Torque
    void setTorque(double torque_){Torque = torque_;}
    
    //! getter for TanglotronMotorUnits
    const std::vector<TanglotronMotor>& getTanglotronMotors() const{return TanglotronMotorUnits;}
    
    /**
    * add one tanglotron motor to TanglotronMotorUnits
    * 
    * @param oneTanglotronMotor_ one TanglotronMotor unit
    * contains indices of one stator and four rotors
    */
    void addTanglotronMotor(TanglotronMotor oneTanglotronMotor_)
    {
        TanglotronMotorUnits.push_back(oneTanglotronMotor_);
    }
    
    //! setter reference for TanglotronMotorUnits
    std::vector<TanglotronMotor>& modifyTanglotronMotors(){return TanglotronMotorUnits;}
    
    /**
    * returns size() of TanglotronMotorUnits
    */
    const size_t getTanglotronMotorsSize() const
    {
        return TanglotronMotorUnits.size();
    }
    
    template<class IngredientsType> 
    void exportRead(FileImport<IngredientsType>& fileReader);
    
    template<class IngredientsType> 
    void exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const;
    
private:
    //! to decide if potential should be calculated or not
    bool variableCalculatePotential;
    
    //! torque
    double Torque;
    //! indices of stators and rotors of all tanglotron motors in ingredients
    std::vector<TanglotronMotor> TanglotronMotorUnits;
    
    template<class IngredientsType> 
    double calculatePotential(const IngredientsType& ingredients, MoveLocalSc& move);
    
    int sgn(double value);
};

/**
* @brief checkMove by probability from calculatePotential
* 
* checks if move is allowed by the Metropolis-criterion.
* 
* @param move MoveLocalSc move
* 
* @tparam IngredientsType Features used in the system. See Ingredients
*/
template<class IngredientsType>
bool FeatureTanglotron::checkMove(const IngredientsType& ingredients, MoveLocalSc& move) 
{	
    if((ingredients.getMolecules()[move.getIndex()].getIsTanglotron()) && (ingredients.getMolecules()[move.getIndex()].getTanglotronType()!=TanglotronAttributeTag::stator)){
        if(ingredients.getCalculatePotential())
            move.multiplyProbability(calculatePotential(ingredients, move));
    }
    
	return true;
}

/**
* @brief calculate potential for Metropolis-criterion
* 
* Look if moved monomer is a rotor monomer and save position of it as well as
* position of stator and the corresponding second rotor monomer.
* Then calculate vectors between stator and each rotor, axis, radius before and
* after move in plane perpendicular to axis and angle between this two radii.
* Then calculate Metropolis-criterion with angle and torque and return it.
* 
* @param move MoveLocalSc move
* 
* @tparam IngredientsType Features used in the system. See Ingredients
*/
template<class IngredientsType> 
double FeatureTanglotron::calculatePotential(const IngredientsType& ingredients, MoveLocalSc& move)
{
    #ifdef DEBUG
        std::cout << "begin to calculate potential" << std::endl;
    #endif /*DEBUG*/
    
    //index of moved monomer
    uint32_t movedMonomerIndex(move.getIndex());
    //ID of moved tanglotron motor
    uint32_t motorID(ingredients.getMolecules()[movedMonomerIndex].getTanglotronID());
    //we need two rotors -> use a vector
    std::vector <VectorInt3> rotorsPosition;
    //we need one stator
    VectorInt3 statorPosition(ingredients.getMolecules()[ingredients.getTanglotronMotors()[motorID].getIndexStator()]);
    
    
    //add first moved rotor monomer to rotorsPosition
    rotorsPosition.push_back(ingredients.getMolecules()[movedMonomerIndex]);
    
    //add the second rotor monomer belongig to the moved rotor monomer
    switch(ingredients.getMolecules()[movedMonomerIndex].getTanglotronType()) // either 1=rotorA, 2=rotorB, 3=rotorC, 4=rotorD
    {
        case TanglotronAttributeTag::rotorA:
            rotorsPosition.push_back(ingredients.getMolecules()[ingredients.getTanglotronMotors()[motorID].getIndexRotorB()]);
            break;
        
        case TanglotronAttributeTag::rotorB:
            rotorsPosition.push_back(ingredients.getMolecules()[ingredients.getTanglotronMotors()[motorID].getIndexRotorA()]);
            break;
        
        case TanglotronAttributeTag::rotorC:
            rotorsPosition.push_back(ingredients.getMolecules()[ingredients.getTanglotronMotors()[motorID].getIndexRotorD()]);
            break;
        
        case TanglotronAttributeTag::rotorD:
            rotorsPosition.push_back(ingredients.getMolecules()[ingredients.getTanglotronMotors()[motorID].getIndexRotorC()]);
            break;
        
        default:
            throw std::runtime_error("FeatureTanglotron: moved monomer is no rotor monomer!");
    }
    //if still here we have two rotor positions in rotorsPosition and one statorPosition and can calculate the potential
    
    //vectors stator->rotor
    //moved vector (according to vector b)
    VectorDouble3 statorRotor1(rotorsPosition.at(0)-statorPosition);
    //non-moved vector (according to vector a)
    VectorDouble3 statorRotor2(rotorsPosition.at(1)-statorPosition);
    
    //axis along statorRotor2 according to normalized vector a
    VectorDouble3 axis(statorRotor2.normalize());
    
    //radius of rotation of vector statorRotor1 according to vector r
    VectorDouble3 radius(statorRotor1-(statorRotor1*axis)*axis);
    //new radius according to vector r' (always constant rotational axis along a)
    VectorDouble3 radiusNew((statorRotor1+move.getDir())-((statorRotor1+move.getDir())*axis)*axis);
    
    //angle between old and new radius
    double angleCos((radius*radiusNew)/(radius.getLength()*radiusNew.getLength()));
    double angle;
    if(radiusNew.getLength()==0) //case for b'||omega => b' comes up to the constant rotationaly axis, but didn't circle in the plane perpendicular to the axis
        return 1.0; //angle = 0;
    else if(angleCos>1.0) //makes sure, that arccos isn't higher than 1 for rounding errors in case of r=r'
            return 1.0; //angle = 0; //acos(1.0);
    else if(angleCos<-1.0) //makes sure, that arccos isn't higher than 1 for rounding errors in case of r=-r'
            angle=M_PI; //acos(-1.0);
    else
        angle=acos(angleCos); //in every non-special case calculate angle according to Delta phi
    
    //normal on plane made of vectors statorRotor2 and statorRotor1 according to vector n=(axb)
    VectorDouble3 normal(crossProduct(statorRotor2,statorRotor1).normalize());
    
    //calculate potential according e^-Delta U
    return (exp(-(sgn(move.getDir()*normal) * Torque * angle)));
}

/**
* calculate signum function
* 
* @param value value to decide if it is above or below 0
*/
int FeatureTanglotron::sgn(double value)
{
    if(value>0)
        return 1;
    else if(value==0)
        return 0;
    else
        return -1;
}

/**
* read value of torque and indices of tanglotron motors from bfm file
* by function ReadTorque and ReadTanglotron
* 
* @tparam IngredientsType Features used in the system. See Ingredients
*/
template < class IngredientsType >
void FeatureTanglotron::exportRead(FileImport<IngredientsType>& fileReader)
{
    IngredientsType& destination=fileReader.getDestination();
    fileReader.registerRead("#!torque", new  ReadTorque< IngredientsType > (destination));
    fileReader.registerRead("#!tanglotrons", new  ReadTanglotron< IngredientsType > (destination));
}

/**
* write value of torque and indices of tanglotron motors in bfm file
* by function WriteTanglotronTorque
* 
* @tparam IngredientsType Features used in the system. See Ingredients
*/
template < class IngredientsType >
void FeatureTanglotron::exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const
{
    const IngredientsType& source=fileWriter.getIngredients_();
    fileWriter.registerWrite("#!torque", new WriteTorque <IngredientsType> (source));
    fileWriter.registerWrite("#!tanglotrons", new WriteTanglotron <IngredientsType> (source));
}

#endif //LEMONADE_FEATURE_TANGLOTRON_H
