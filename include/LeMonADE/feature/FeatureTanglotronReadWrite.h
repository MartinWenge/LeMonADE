 /*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2021 by 
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

#ifndef FEATURE_TANGLOTRON_READ_WRITE_H
#define FEATURE_TANGLOTRON_READ_WRITE_H
/**
 * @file
 * @date 2021/02/02
 * @author Martin Wengenmayr/Cornelia Struebig
 * @brief Definition and Implementation of Read/Write utility of FeatureTanglotron
**/

#include <iostream>
#include <sstream>

#include <LeMonADE/io/AbstractRead.h>

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
            for(uint n=0; n<(this->getSource()).getTanglotronMotors().size(); n++){
                strm << (this->getSource()).getTanglotronMotors()[n].getIndexStator()+1 << " ";
                strm << (this->getSource()).getTanglotronMotors()[n].getIndexRotorA()+1 << " ";
                strm << (this->getSource()).getTanglotronMotors()[n].getIndexRotorB()+1 << " ";
                strm << (this->getSource()).getTanglotronMotors()[n].getIndexRotorC()+1 << " ";
                strm << (this->getSource()).getTanglotronMotors()[n].getIndexRotorD()+1 << "\n";
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
            TanglotronMotor newTanglotron;
            newTanglotron.setTanglotronMotor(idxStat-1, idxA-1, idxB-1, idxC-1, idxD-1);
            this->getDestination().addTanglotronMotor(newTanglotron);
            
            this->getDestination().modifyMolecules()[idxStat-1].setIsTanglotron(true);
            this->getDestination().modifyMolecules()[idxStat-1].setTanglotronType(TanglotronAttributeTag::tanglotronType::stator);
            this->getDestination().modifyMolecules()[idxStat-1].setTanglotronID(this->getDestination().getTanglotronMotorsSize()-1);
            this->getDestination().modifyMolecules()[idxA-1].setIsTanglotron(true);
            this->getDestination().modifyMolecules()[idxA-1].setTanglotronType(TanglotronAttributeTag::tanglotronType::rotorA);
            this->getDestination().modifyMolecules()[idxA-1].setTanglotronID(this->getDestination().getTanglotronMotorsSize()-1);
            this->getDestination().modifyMolecules()[idxB-1].setIsTanglotron(true);
            this->getDestination().modifyMolecules()[idxB-1].setTanglotronType(TanglotronAttributeTag::tanglotronType::rotorB);
            this->getDestination().modifyMolecules()[idxB-1].setTanglotronID(this->getDestination().getTanglotronMotorsSize()-1);
            this->getDestination().modifyMolecules()[idxC-1].setIsTanglotron(true);
            this->getDestination().modifyMolecules()[idxC-1].setTanglotronType(TanglotronAttributeTag::tanglotronType::rotorC);
            this->getDestination().modifyMolecules()[idxC-1].setTanglotronID(this->getDestination().getTanglotronMotorsSize()-1);
            this->getDestination().modifyMolecules()[idxD-1].setIsTanglotron(true);
            this->getDestination().modifyMolecules()[idxD-1].setTanglotronType(TanglotronAttributeTag::tanglotronType::rotorD);
            this->getDestination().modifyMolecules()[idxD-1].setTanglotronID(this->getDestination().getTanglotronMotorsSize()-1);
            
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


#endif //FEATURE_TANGLOTRON_READ_WRITE_H
