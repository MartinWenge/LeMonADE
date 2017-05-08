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


#ifndef LEMONADE_FEATUREADSORPTION_H
#define LEMONADE_FEATUREADSORPTION_H

/**
 * @file
 * @date 2017/05/05
 * @author Martin Wengenmayr
 * @brief Definition and implementation of class FeatureAdsorption
**/

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/feature/FeatureBoltzmann.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/feature/FeatureAdsorptionReadWrite.h>

#include <cmath> // for std::exp

/**
 * @class FeatureAdsorption
 * @brief implementation of adsorption 
 *
 * @tparam IngredientsType type of ingredients
 *
 * @details
 * Adsorption energie per monomer is used to calculate the move probability if monoer is on the adsorbing wall
**/

class FeatureAdsorption:public Feature
{

private:
  //! boolean for x-Wall adsorption
  bool adsorptionX;
  
  //! boolean for y-Wall adsorption
  bool adsorptionY;
  
  //! boolean for z-Wall adsorption
  bool adsorptionZ;
  
  //! adsorption energy for x-Wall adsorption
  double adsorptionEnergyX;
  
  //! adsorption energy for y-Wall adsorption
  double adsorptionEnergyY;
  
  //! adsorption energy for z-Wall adsorption
  double adsorptionEnergyZ;
  
  //! probability calculated from Boltzmann criteria and adsorption energy for x-Wall adsorption
  double probabilityX;
  
  //! probability calculated from Boltzmann criteria and adsorption energy for y-Wall adsorption
  double probabilityY;
  
  //! probability calculated from Boltzmann criteria and adsorption energy for z-Wall adsorption
  double probabilityZ;
  
  //! checker if adsorption energy is set
  bool energyIsSet;
 

public:

  FeatureAdsorption();
  ~FeatureAdsorption(){}
    
  //! required back feature: FeatureBoltzmann for applying Metropolis criterion
  typedef LOKI_TYPELIST_1(FeatureBoltzmann) required_features_back;

  //! required front feature: FeatureBox
  typedef LOKI_TYPELIST_1(FeatureBox) required_features_front;

  //! check for all Monte Carlo moves without special check functions (always true)
  template<class IngredientsType>
  bool checkMove(const IngredientsType& ingredients, const MoveBase& move) const;

  //! check for local moves
  template<class IngredientsType,class LocalMoveType> 
  bool checkMove(const IngredientsType& ingredients, MoveLocalBase<LocalMoveType>& move);
  
  //! synchronize function to check if Box settings match the adsorption settings
  template < class IngredientsType > void synchronize(IngredientsType& ing);

  /* ********* X direction ******** */
  //! getter function of adsorption energie in x direction
  double getAdsorptionEnergyX() const{return adsorptionEnergyX;}
  //! setter for adsorption energie in x direction
  void setAdsorptionEnergyX(double ae_){
    adsorptionEnergyX=ae_;
    energyIsSet=true;
    probabilityX=std::exp(-ae_);
    adsorptionX=true;
  }
  //! getter function of adsorption bool in x direction
  bool getAdsorptionX() const{return adsorptionX;}
  //! setter function of adsorption bool in x direction
  void setAdsorptionX(bool adsorpX_){adsorptionX=adsorpX_;}
  
  /* ********* Y direction ******** */
  //! getter function of adsorption energie in y direction
  double getAdsorptionEnergyY() const{return adsorptionEnergyY;}
  //! setter for adsorption energie in y direction
  void setAdsorptionEnergyY(double ae_){
    adsorptionEnergyY=ae_;
    energyIsSet=true;
    probabilityY=std::exp(-ae_);
    adsorptionY=true;
  }
  //! getter function of adsorption bool in y direction
  bool getAdsorptionY() const{return adsorptionY;}
  //! setter function of adsorption bool in y direction
  void setAdsorptionY(bool adsorpY_){adsorptionY=adsorpY_;}
  
  /* ********* Z direction ******** */
  //! getter function of adsorption energie in z direction
  double getAdsorptionEnergyZ() const{return adsorptionEnergyZ;}
  //! setter for adsorption energie in z direction
  void setAdsorptionEnergyZ(double ae_){
    adsorptionEnergyZ=ae_;
    energyIsSet=true;
    probabilityZ=std::exp(-ae_);
    adsorptionZ=true;
  }
  //! getter function of adsorption bool in z direction
  bool getAdsorptionZ() const{return adsorptionZ;}
  //! setter function of adsorption bool in z direction
  void setAdsorptionZ(bool adsorpZ_){adsorptionZ=adsorpZ_;}
  

  //!export bfm-file read command !nn_interaction
  template <class IngredientsType>
  void exportRead(FileImport <IngredientsType>& fileReader);

  //!export bfm-file write command !nn_interaction
  template <class IngredientsType>
  void exportWrite(AnalyzerWriteBfmFile <IngredientsType>& fileWriter) const;
  
  void printMetaData(std::ostream& streamOut) const{
    streamOut << "\tAbdsorption in (x,y,z):("  << adsorptionX<<","<< adsorptionY<<"," << adsorptionZ<<")"<< std::endl;
  }

};


/**
 * @brief Constructor
 **/
FeatureAdsorption::FeatureAdsorption():
  adsorptionX(false), adsorptionY(false), adsorptionZ(false),
  adsorptionEnergyX(0.0), adsorptionEnergyY(0.0),adsorptionEnergyZ(0.0),
  probabilityX(1.0), probabilityY(1.0), probabilityZ(1.0), 
  energyIsSet(false)
{}

/**
 * @details Returns true for all moves other than the ones that have specialized versions of this function.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move Monte Carlo move other than moves with specialized functions.
 * @return true (always)
 **/
template<class IngredientsType>
bool FeatureAdsorption::checkMove(const IngredientsType& ingredients, const MoveBase& move) const {
    return true;
}

/**
 * @details calculates the factor for the acceptance probability of the move
 * arising from the monomer position at the wall or away from it
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move Monte Carlo move of type MoveLocalSc
 * @return true (always)
 * 
 * @details
 * Find all possible adsorption contacts, and add the probability factor 
 * coming from this feature, then return true,
 * because the total probability is evaluated by FeatureBoltzmann at the end.
 * 
 **/
 
template<class IngredientsType,class LocalMoveType> 
bool FeatureAdsorption::checkMove(const IngredientsType& ingredients, MoveLocalBase<LocalMoveType>& move) {
  VectorInt3 dir(move.getDir());
  
  // check adsorption on X-Walls 
  if( adsorptionX && (dir.getX()!=0) ){
    int32_t pos(ingredients.getMolecules()[move.getIndex()].getX());
    if( (pos==0 && dir.getX()==1) || 
      (pos==(ingredients.getBoxX()-1) && dir.getX()==-1) 
    ){
      // decrease move probability if monomer whants to leave the wall
      move.multiplyProbability(probabilityX);
    }else if( (pos==1 && dir.getX()==-1) || 
      (pos==(ingredients.getBoxX()-2) && dir.getX()==1) 
    ){
      // increase the move probability if monomer wants to touch the wall
      move.multiplyProbability(1.0/probabilityX);
    }
  }
    
  // check adsorption on Y-Walls 
  if( adsorptionY && (dir.getY()!=0) ){
    int32_t pos(ingredients.getMolecules()[move.getIndex()].getY());
    if( (pos==0 && dir.getY()==1) || 
      (pos==(ingredients.getBoxY()-1) && dir.getY()==-1) 
    ){
      // decrease move probability if monomer whants to leave the wall
      move.multiplyProbability(probabilityY);
    }else if( (pos==1 && dir.getY()==-1) || 
      (pos==(ingredients.getBoxY()-2) && dir.getY()==1) 
    ){
      // increase the move probability if monomer wants to touch the wall
      move.multiplyProbability(1.0/probabilityY);
    }
  }
  
  // check adsorption on X-Walls 
  if( adsorptionZ && (dir.getZ()!=0) ){
    int32_t pos(ingredients.getMolecules()[move.getIndex()].getZ());
    if( (pos==0 && dir.getZ()==1) || 
      (pos==(ingredients.getBoxZ()-1) && dir.getZ()==-1) 
    ){
      // decrease move probability if monomer whants to leave the wall
      move.multiplyProbability(probabilityZ);
    }else if( (pos==1 && dir.getZ()==-1) || 
      (pos==(ingredients.getBoxZ()-2) && dir.getZ()==1) 
    ){
      // increase the move probability if monomer wants to touch the wall
      move.multiplyProbability(1.0/probabilityZ);
    }
  }
  
  return true;
}

/**
 * @tparam IngredientsType The type of the system including all features
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * 
 * @details Check if FeatureWall setting match FeatureAdsorption settings
 * and if adsorptionEnergy was set.
 **/
 
template<class IngredientsType>
void FeatureAdsorption::synchronize(IngredientsType& ingredients){
  if(adsorptionX){
    if(ingredients.isPeriodicX())
      throw std::runtime_error("FeatureAdsorption: adsorption on periodic X walls not allowed!");
  }
  if(adsorptionY){
    if(ingredients.isPeriodicY())
      throw std::runtime_error("FeatureAdsorption: adsorption on periodic Y walls not allowed!");
  }
  if(adsorptionZ){
    if(ingredients.isPeriodicZ())
      throw std::runtime_error("FeatureAdsorption: adsorption on periodic Z walls not allowed!");
  }
  // check if one of the energies is set 
  //! @todo too restrictive?
  if(!energyIsSet)
    throw std::runtime_error("FeatureAdsorption: adsorption energy is not set!");
}


/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type FileImport. The export of the Reads is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * - !adsorption
 *
 * @tparam IngredientsType The type of the system including all features
 * @param fileReader File importer for the bfm-file
 **/
 
template<class IngredientsType>
void FeatureAdsorption::exportRead(FileImport< IngredientsType >& fileReader)
{
  fileReader.registerRead("!adsorption",new ReadAdsorption< FeatureAdsorption >(*this));
}

/**
 * The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Write-Out Commands:
 * - !adsorption
 *
 * @tparam IngredientsType The type of the system including all features
 * @param fileWriter File writer for the bfm-file.
 */
 
template<class IngredientsType>
void FeatureAdsorption::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& fileWriter) const
{
  fileWriter.registerWrite("!adsorption",new WriteAdsorption< FeatureAdsorption >(*this));
}


#endif /*LEMONADE_FEATUREADSORPTION_H*/
