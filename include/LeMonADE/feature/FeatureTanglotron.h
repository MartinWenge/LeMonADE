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
#include <map>

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/io/AbstractRead.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/feature/FeatureBoltzmann.h>

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
   
	//! Standard constructor- initially the tag is set to -1 and isTanglotron is set to false
	TanglotronAttributeTag():isTanglotron(false),tag(-1),tanglotronID(-1){}

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
	int32_t getTanglotronType() const {return tag;}
	
	/**
	* @brief setting size_t tag tanglotronID
	*
	* @param id number of motor in vector TanglotronMotorUnits
	**/
    void setTanglotronID(int32_t id){ tanglotronID=id; }
    
    //! getting size_t tag tanglotronID
	int32_t getTanglotronID() const {return tanglotronID;}
	
private:
    //! Private Variable holding the bool if monomer is part of tanglotron
    bool isTanglotron;
	//! Private variable holding the tanglotronType. Default is 5.
    int32_t tag;
    //! private variable holding the number of one tanglotron motor in the vector holding all motors
    int32_t tanglotronID;
};


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
    FeatureTanglotron():variableCalculatePotential(true), Torque(0.0) {}
    
    //! default empty destructor
    virtual ~FeatureTanglotron(){}
    
    /**
    * synchronize
    *
    * @param ingredients A reference to the IngredientsType - mainly the system.
    */
    template<class IngredientsType> void synchronize(IngredientsType& ingredients)
    {
        if(tanglotronMotorMap.empty()){
            tanglotronMotorMap.clear();
        }
        // loop over the system to find tanglotrons
        for(uint32_t i=0; i<ingredients.getMolecules().size(); i++){
            // check if monomer is part of a tanglotron
            if(ingredients.getMolecules()[i].getIsTanglotron()){
                // check if tanglotronID was already added to the tanglotronMotorMap
                int32_t motorID(ingredients.getMolecules()[i].getTanglotronID());
                auto map_it = tanglotronMotorMap.find(motorID);
                // if motorID is new, add new vector to the map
                if(map_it == tanglotronMotorMap.end()){
                    tanglotronMotorMap[motorID] = std::vector<int32_t>(5,0);
                }
                // add the tanglotron type to the vector
                tanglotronMotorMap[motorID].at(ingredients.getMolecules()[i].getTanglotronType()) = i;
            }
        }
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

    //! getter for TanglotronMap
    std::map<int32_t, std::vector <int32_t> > getTanglotronMotorMap() const{return tanglotronMotorMap;}
    
    template<class IngredientsType> 
    void exportRead(FileImport<IngredientsType>& fileReader);
    
    template<class IngredientsType> 
    void exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const;
    
private:
    //! to decide if potential should be calculated or not
    bool variableCalculatePotential;
    
    //! torque
    double Torque;

    //! lookup for tanglotron motors
    std::map<int, std::vector<int32_t> > tanglotronMotorMap;
    
    template<class IngredientsType> 
    double calculatePotential(const IngredientsType& ingredients, MoveLocalSc& move);
    
    template <typename T> int sgn(T val) {return (T(0) < val) - (val < T(0));}
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
    VectorInt3 statorPosition(ingredients.getMolecules()[tanglotronMotorMap.at(motorID).at(TanglotronAttributeTag::tanglotronType::stator)]);
    
    //add first moved rotor monomer to rotorsPosition
    rotorsPosition.push_back(ingredients.getMolecules()[movedMonomerIndex]);
    
    //add the second rotor monomer belonging to the moved rotor monomer
    switch(ingredients.getMolecules()[movedMonomerIndex].getTanglotronType()) // either 1=rotorA, 2=rotorB, 3=rotorC, 4=rotorD
    {
        case TanglotronAttributeTag::rotorA:
            rotorsPosition.push_back(ingredients.getMolecules()[tanglotronMotorMap.at(motorID).at(TanglotronAttributeTag::tanglotronType::rotorB)]);
            break;
        
        case TanglotronAttributeTag::rotorB:
            rotorsPosition.push_back(ingredients.getMolecules()[tanglotronMotorMap.at(motorID).at(TanglotronAttributeTag::tanglotronType::rotorA)]);
            break;
        
        case TanglotronAttributeTag::rotorC:
            rotorsPosition.push_back(ingredients.getMolecules()[tanglotronMotorMap.at(motorID).at(TanglotronAttributeTag::tanglotronType::rotorD)]);
            break;
        
        case TanglotronAttributeTag::rotorD:
            rotorsPosition.push_back(ingredients.getMolecules()[tanglotronMotorMap.at(motorID).at(TanglotronAttributeTag::tanglotronType::rotorC)]);
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
    if(radiusNew.getLength()==0) //case for b'||omega => b' comes up to the constant rotationally axis, but didn't circle in the plane perpendicular to the axis
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
*
* @class ReadTanglotron
*
* @brief reading #!tanglotrons from bfm file
* 
* @details Form:
* #!tanglotrons
* indexStator indexRotorA indexRotorB indexRotorC indexRotorD
*
* @tparam IngredientsType
**/

template<class IngredientsType>
class ReadTanglotron: public ReadToDestination<IngredientsType>
{
    public:
        //! constructor
        ReadTanglotron(IngredientsType& destination): ReadToDestination<IngredientsType>(destination){}

        void execute();
};

/**
* @class WriteTanglotron
*
* @brief writing #!tanglotrons to bfm file
* 
* @details Form:
* #!tanglotrons
* indexStator indexRotorA indexRotorB indexRotorC indexRotorD
*
* @tparam IngredientsType
**/
template<class IngredientsType>
class WriteTanglotron: public AbstractWrite<IngredientsType>
{
    public:
        //! constructor
        WriteTanglotron(const IngredientsType& src):AbstractWrite<IngredientsType>(src){this->setHeaderOnly(true);}
  
        /**
        * writeStream
        * 
        * @brief getTanglotronMotors and write it to bfm file
        */ 
        void writeStream(std::ostream& strm)
        {
            const IngredientsType& ingredients=(this->getSource());
            
            //add 1 to indices of tanglotron monomers because of different index
            //definitions of bfm and LeMonADE
            strm << "#!tanglotrons\n";

            std::map<int, std::vector<int32_t> > tanglotronMotorMap;
            // loop over the system to find tanglotrons
            for(uint32_t i=0; i<ingredients.getMolecules().size(); i++){
                // check if monomer is part of a tanglotron
                if(ingredients.getMolecules()[i].getIsTanglotron()){
                    // check if tanglotronID was already added to the tanglotronMotorMap
                    int32_t motorID(ingredients.getMolecules()[i].getTanglotronID());
                    auto map_it = tanglotronMotorMap.find(motorID);
                    // if motorID is new, add new vector to the map
                    if(map_it == tanglotronMotorMap.end()){
                        tanglotronMotorMap[motorID] = std::vector<int32_t>(5,0);
                    }
                    // add the tanglotron type to the vector
                    tanglotronMotorMap[motorID].at(ingredients.getMolecules()[i].getTanglotronType()) = i;
                }
            }

            for(auto const& motor: tanglotronMotorMap){
                strm << motor.second.at(0)+1 << " ";
                strm << motor.second.at(1)+1 << " ";
                strm << motor.second.at(2)+1 << " ";
                strm << motor.second.at(3)+1 << " ";
                strm << motor.second.at(4)+1 << "\n";
            }
            strm << "\n";
        }
};

/**
* @class ReadTorque
*
* @brief reading #!torque from bfm file
* 
* @details Form: #!torque=10
* 
* @tparam IngredientsType
**/
template<class IngredientsType>
class ReadTorque: public ReadToDestination<IngredientsType>
{
    public:
        //! constructor
        ReadTorque(IngredientsType& destination):ReadToDestination<IngredientsType>(destination){}
        
        void execute();
};

/**
* @class WriteTorque
*
* @brief writing #!torque to bfm file
* 
* @details Form: #!torque=10
*
* @tparam IngredientsType
**/
template<class IngredientsType>
class WriteTorque: public AbstractWrite<IngredientsType>
{
    public:
        //! constructor
        WriteTorque(const IngredientsType& src):AbstractWrite<IngredientsType>(src){this->setHeaderOnly(true);}
  
        /**
        * writeStream
        * 
        * @brief getTorque and write it to bfm file
        */ 
        void writeStream(std::ostream& strm)
        {
            const IngredientsType& ingredients=(this->getSource());
            
            strm << "#!torque=" << ingredients.getTorque() << "\n\n";
        }
};


/////////////MEMBER IMPLEMENTATIONS ////////////////////////////////////////////
/**
* ReadTanglotron::execute()
* 
* @brief reading tanglotrons and add them to private vector TanglotronMotorUnits
*/
template<class IngredientsType>
void ReadTanglotron<IngredientsType>::execute()
{
    std::cout << "reading #!tanglotrons ...\n";
    
    uint32_t idxStat, idxA, idxB, idxC, idxD;
    int32_t tanglotronID(0);
    std::string line;
    std::streampos previous;
    
    std::getline(this->getInputStream(),line);
    previous=(this->getInputStream()).tellg();
    getline(this->getInputStream(),line);
    
    while(!line.empty() && !((this->getInputStream()).fail())){

        //stop at next Read and set the get-pointer to the position before the Read
        if(this->detectRead(line)){
            (this->getInputStream()).seekg(previous);
            break;
        }
  
        std::stringstream stream(line);
        
        //read index of stator, throw exception if extraction fails
        stream >> idxStat >> idxA >> idxB >> idxC >> idxD;
        if(stream.fail()){
            std::stringstream messagestream;
            messagestream << "ReadTanglotron<IngredientsType>::execute() -> Could not read indices";
            throw std::runtime_error(messagestream.str());
        }
        

        //if streaming worked set the TanglotronMotor with indices -1 because of
        //different index definitions of bfm and LeMonADE and add ist to private
        //vector TanglotronMotorUnits
        if(!stream.fail()){
            this->getDestination().modifyMolecules()[idxStat-1].setIsTanglotron(true);
            this->getDestination().modifyMolecules()[idxStat-1].setTanglotronType(TanglotronAttributeTag::tanglotronType::stator);
            this->getDestination().modifyMolecules()[idxStat-1].setTanglotronID(tanglotronID);
            this->getDestination().modifyMolecules()[idxA-1].setIsTanglotron(true);
            this->getDestination().modifyMolecules()[idxA-1].setTanglotronType(TanglotronAttributeTag::tanglotronType::rotorA);
            this->getDestination().modifyMolecules()[idxA-1].setTanglotronID(tanglotronID);
            this->getDestination().modifyMolecules()[idxB-1].setIsTanglotron(true);
            this->getDestination().modifyMolecules()[idxB-1].setTanglotronType(TanglotronAttributeTag::tanglotronType::rotorB);
            this->getDestination().modifyMolecules()[idxB-1].setTanglotronID(tanglotronID);
            this->getDestination().modifyMolecules()[idxC-1].setIsTanglotron(true);
            this->getDestination().modifyMolecules()[idxC-1].setTanglotronType(TanglotronAttributeTag::tanglotronType::rotorC);
            this->getDestination().modifyMolecules()[idxC-1].setTanglotronID(tanglotronID);
            this->getDestination().modifyMolecules()[idxD-1].setIsTanglotron(true);
            this->getDestination().modifyMolecules()[idxD-1].setTanglotronType(TanglotronAttributeTag::tanglotronType::rotorD);
            this->getDestination().modifyMolecules()[idxD-1].setTanglotronID(tanglotronID);

            tanglotronID++;
            
            std::getline(this->getInputStream(),line);
        }
        
        //otherwise throw an exception
        else{
            std::stringstream messagestream;
            messagestream << "ReadTanglotron<IngredientsType>::execute()\n" << "Could not read indices in #!tanglotrons";
            throw std::runtime_error(messagestream.str());
    
        }
    }    
}

/**
* ReadTorque::execute()
* 
* @brief reading torque and set private variable Torque
*/
template<class IngredientsType>
void ReadTorque<IngredientsType>::execute()
{
    std::cout << "reading #!torque...";
    
    double torque;
    
    if(this->getInputStream() >> torque)
    { 
        this->getDestination().setTorque(torque);
        std::cout << torque << std::endl;
    }
    else
        throw std::runtime_error("ReadTorque<IngredientsType>::execute()\n Could not read torque");    
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
