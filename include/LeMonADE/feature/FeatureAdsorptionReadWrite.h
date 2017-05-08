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

#ifndef LEMONADE_FEATUREADSORPTION_READ_WRITE_H
#define LEMONADE_FEATUREADSORPTION_READ_WRITE_H

/**
 * @file
 * @date 2017/05/05
 * @author Martin wengenmayr
 * @brief Def. and impl. of class templates ReadAdsorption and WriteAdsorption
**/

#include<iostream>
#include<LeMonADE/io/AbstractRead.h>
#include<LeMonADE/io/AbstractWrite.h>

/**
 * @class ReadAdsorption
 * @brief Handles BFM-file read command !nn_interaction
 * @tparam IngredientsType Ingredients class storing all system information.
**/
template < class IngredientsType>
class ReadAdsorption: public ReadToDestination<IngredientsType>
{
public:
    ReadAdsorption(IngredientsType& ingredients):
      ReadToDestination<IngredientsType>(ingredients){}
    virtual ~ReadAdsorption(){}
    virtual void execute();
};

/**
 * @class WriteAdsorption
 * @brief Handles BFM-file write command !nn_interaction
 * @tparam IngredientsType Ingredients class storing all system information.
**/
template <class IngredientsType>
class WriteAdsorption:public AbstractWrite<IngredientsType>
{
public:

  //constructor sets the headerOnly tag such that
  // adsorption information is only written once in the header
    WriteAdsorption(const IngredientsType& ingredients):
      AbstractWrite<IngredientsType>(ingredients){this->setHeaderOnly(true);}

    virtual ~WriteAdsorption(){}

    virtual void writeStream(std::ostream& strm);
};

/////////////MEMBER IMPLEMENTATIONS ////////////////////////////////////////////

/**
 * @brief Executes the reading routine to extract \b !nn_interaction.
 *
 * @throw <std::runtime_error> fail to read monomer types or interaction energy.
 **/
template<class IngredientsType>
void ReadAdsorption<IngredientsType>::execute()
{
  
    IngredientsType& ingredients=this->getDestination();
    std::istream& file=this->getInputStream();

    double adsorpEnergy;
    std::string direction;

    //set stream to throw exception on fail
    file.exceptions(file.exceptions() | std::ifstream::failbit);
    
    try
      {
	file>>adsorpEnergy;
	file>>direction;
      }
    catch(std::ifstream::failure e)
      {
	std::stringstream errormessage;
	errormessage<<"ReadAdsorption::execute().\n";
	errormessage<<"Could not read adsorption from file\n";
	errormessage<<"Previous error: "<<e.what()<<std::endl;
	throw std::runtime_error(errormessage.str());
      }
    
    //now save the interaction tuple just read from the file
    if( direction.compare("x")==0 ){
      ingredients.setAdsorptionEnergyX(adsorpEnergy);
    }else if( direction.compare("y")==0 ){
      ingredients.setAdsorptionEnergyY(adsorpEnergy);
    }else if( direction.compare("z")==0 ){
      ingredients.setAdsorptionEnergyZ(adsorpEnergy);
    }else{
      throw std::runtime_error("ReadAdsorption::execute(): unknown read in key\n");
    }

    //unset exception on fail
    file.exceptions(file.exceptions() & (~std::ifstream::failbit));
}


/**
 * @brief Executes the routine to write \b !nn_interaction.
 * @arg stream file stream to write into
 **/
template<class IngredientsType>
void WriteAdsorption<IngredientsType>::writeStream(std::ostream& stream)
{
  stream<<"## adsorbing walls with respective adsorption energy\n";
  if(this->getSource().getAdsorptionX()){
    stream<<"!adsorption "<<this->getSource().getAdsorptionEnergyX()<<" x\n";
  }
  if(this->getSource().getAdsorptionY()){
    stream<<"!adsorption "<<this->getSource().getAdsorptionEnergyY()<<" y\n";
  }
  if(this->getSource().getAdsorptionZ()){
    stream<<"!adsorption "<<this->getSource().getAdsorptionEnergyZ()<<" z\n";
  }
  
  stream<<"\n\n";

}




#endif // LEMONADE_FEATUREADSORPTION_READ_WRITE_H
