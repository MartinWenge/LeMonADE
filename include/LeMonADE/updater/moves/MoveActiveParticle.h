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

#ifndef LEMONADE_UPDATER_MOVES_MOVELOCALACTIVE_H
#define LEMONADE_UPDATER_MOVES_MOVELOCALACTIVE_H

#include <LeMonADE/updater/moves/MoveLocalBase.h>

/*****************************************************************************/
/**
 * @file
 *
 * @class MoveActiveParticle
 *
 * @brief move modeling "active" particle
 *
 * @details The class is a specialization of MoveLocalBase using the (CRTP) to avoid virtual functions.
 **/
/*****************************************************************************/
class MoveActiveParticle:public MoveLocalBase<MoveActiveParticle>
{
public:
  MoveActiveParticle(){
	steps[0]=VectorInt3(1,0,0);
	steps[1]=VectorInt3(-1,0,0);
	steps[2]=VectorInt3(0,1,0);
	steps[3]=VectorInt3(0,-1,0);
	steps[4]=VectorInt3(0,0,1);
	steps[5]=VectorInt3(0,0,-1);	
  }
  
  // overload initialise function to be able to set the moves index and direction if neccessary
  template <class IngredientsType> void init(const IngredientsType& ing);
  template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index);
  
  template <class IngredientsType> bool check(IngredientsType& ing);
  template< class IngredientsType> void apply(IngredientsType& ing);
    
private:
  // holds the possible move directions
  /**
   * @brief Array that holds the 6 possible move directions
   *
   * @details In the scBFM the classic moves (dx,dy,dz) are along the lattice-axes as:
   * * steps   = (dx, dy, dz)
   * * steps[0]= ( 1,  0,  0);
   * * steps[1]= (-1,  0,  0);
   * * steps[2]= ( 0,  1,  0);
   * * steps[3]= ( 0, -1,  0);
   * * steps[4]= ( 0,  0,  1);
   * * steps[5]= ( 0,  0, -1);
   */
  VectorInt3 steps[6];

  template <class IngredientsType>uint32_t getOrientation(const IngredientsType& ing);
};



/////////////////////////////////////////////////////////////////////////////
/////////// implementation of the members ///////////////////////////////////

/*****************************************************************************/
/**
 * @brief Initialize the move.
 *
 * @details Resets the move probability to unity. Dice a new random direction and
 * Vertex (monomer) index inside the graph.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 **/
template <class IngredientsType>
void MoveActiveParticle::init(const IngredientsType& ing)
{
  this->resetProbability();
  
  //draw index
  this->setIndex( (this->randomNumbers.r250_rand32()) %(ing.getMolecules().size()) );
  
  //draw direction
  uint32_t activeDir=getOrientation(ing);
  this->setDir(steps[activeDir]);

}

/*****************************************************************************/
/**
 * @brief Initialize the move with a given monomer index.
 *
 * @details Resets the move probability to unity. Dice a new random direction.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @param index index of the monomer to be moved
 **/
template <class IngredientsType>
void MoveActiveParticle::init(const IngredientsType& ing, uint32_t index)
{
  this->resetProbability();
  
  //set index
  if( (index >= 0) && (index <= (ing.getMolecules().size()-1)) )
    this->setIndex( index );
  else
    throw std::runtime_error("MoveActiveParticle::init(ing, index): index out of range!");
  
  //draw direction
  uint32_t activeDir=getOrientation(ing);
  this->setDir(steps[activeDir]);

}


/*****************************************************************************/
/**
 * @brief Check if the move is accepted by the system.
 *
 * @details This function delegates the checking to the Feature.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @return True if move is valid. False, otherwise.
 **/
template <class IngredientsType>
bool MoveActiveParticle::check(IngredientsType& ing)
{
  //send the move to the Features to be checked
  return ing.checkMove(ing,*this);
}
  
/*****************************************************************************/
/**
 * @brief Apply the move to the system , e.g. add the displacement to Vertex (monomer) position.
 *
 * @details As first step: all Feature should apply the move using applyMove().\n
 * Second: Modify the positions etc. of the Vertex etc.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 **/
template< class IngredientsType>
void MoveActiveParticle::apply(IngredientsType& ing)
{
	///@todo Think about the applying of move. Esp. make this independent of the order to avoid confusion!!
	
	//move must FIRST be applied to the features
	ing.applyMove(ing,*this);	

	//THEN the position can be modified
	ing.modifyMolecules()[this->getIndex()]+=this->getDir();
  
}

/*****************************************************************************/
/**
 * @brief helper function to get average orientation
 *
 * @details average over neighbors and find the ones to orient with
 *
 * @param ing A reference to the IngredientsType - mainly the system
 **/
template <class IngredientsType>
uint32_t MoveActiveParticle::getOrientation(const IngredientsType& ingredients){
// Index of moved Monomer is monoIndex
	uint32_t monoIndex=this->getIndex();
	const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();
	std::map<double,std::vector<uint32_t> > distance_map;

	//quick and dirty: get the n nearest neigbors
	for(uint32_t i=0; i<ingredients.getMolecules().size();i++){
		if(i != monoIndex){
			// get total distance ( maybe replace by squared distance )
      double dx = ingredients.getMolecules()[monoIndex].getX()-ingredients.getMolecules()[i].getX();
      dx= dx - ingredients.getBoxX() * round( dx / ingredients.getBoxX() );
      double dy = ingredients.getMolecules()[monoIndex].getY()-ingredients.getMolecules()[i].getY();
      dy= dy - ingredients.getBoxY() * round( dy / ingredients.getBoxY() );
      double dz = ingredients.getMolecules()[monoIndex].getZ()-ingredients.getMolecules()[i].getZ();
      dz= dz - ingredients.getBoxZ() * round( dz / ingredients.getBoxZ() );
			double distance( dx*dx + dy*dy + dz*dz );
      std::cout << "sort entries: dist= "<< distance<<", idx= "<< i <<std::endl; 
			// if map size is sufficient first check if distance is already larger than the distance of the last particle to consider, else just add it
			if(distance_map.size() > ingredients.getNumOfNeighbors()){
				std::map<double,std::vector<uint32_t> >::iterator it=distance_map.begin();
				std::advance(it,ingredients.getNumOfNeighbors());
				if(distance < ( it->first ) ){
          std::map<double,std::vector<uint32_t> >::iterator findme = distance_map.find(distance);
          if( findme != distance_map.end() ){
            findme->second.push_back(i);
          }else{
            distance_map[distance]=std::vector<uint32_t>(1,i);
          }
				}
			}else{
				std::map<double,std::vector<uint32_t> >::iterator findme = distance_map.find(distance);
        if( findme != distance_map.end() ){
          findme->second.push_back(i);
        }else{
          distance_map[distance]=std::vector<uint32_t>(1,i);
        }
			}
		}
	}

  for( const auto& myMap : distance_map){
    std::cout << myMap.first << "\t";
    for(uint32_t i=0; i< myMap.second.size(); i++){
      std::cout<< myMap.second.at(i)<<" ";
    }
    std::cout << std::endl;
  }

	// get orientation of the nearest neighbors
	VectorDouble3 movesOrientation(ingredients.getMolecules()[monoIndex].getOrientation());
  std::cout << "my   idx = " << monoIndex << " distance = 0"<<std::endl;
	std::map<double,std::vector<uint32_t>>::iterator it=distance_map.begin();
	for( uint32_t i = 0; i < ingredients.getNumOfNeighbors(); i++){
    for( uint32_t j = 0; j < it->second.size(); j++){
      movesOrientation+=ingredients.getMolecules()[it->second.at(j)].getOrientation();
      std::cout << "mono idx = " << it->second.at(j) << " distance = "<< it->first<<std::endl;
      if(j>0){
        i++;
      }
    } 
		it++;
	}
	movesOrientation/=(ingredients.getNumOfNeighbors()+1);

	// add some noise
	VectorDouble3 movesNoise(randomNumbers.r250_drand(),randomNumbers.r250_drand(),randomNumbers.r250_drand());
	movesOrientation+=(ingredients.getNoise() * movesNoise);

	// get the new direction
	double absX=std::abs(movesOrientation.getX());
	double absY=std::abs(movesOrientation.getY());
	double absZ=std::abs(movesOrientation.getZ());

  int32_t dirIdx;
	if( (absX > absY) && (absX > absZ ) ){
		std::cout<<"case 1: "<< absX << " " << absY << " " << absZ << std::endl;
    if(absX > 0){
      dirIdx=0;
    }else{
      dirIdx=1;
    }
	}else if( (absY > absX) && (absY > absZ ) ){
		std::cout<<"case 2: "<< absX << " " << absY << " " << absZ << std::endl;
    if(absX > 0){
      dirIdx=2;
    }else{
      dirIdx=3;
    }
	}else{ // if( (absZ > absX) && (absZ > absY ) ){
		std::cout<<"case 3: "<< absX << " " << absY << " " << absZ << std::endl;
    if(absX > 0){
      dirIdx=4;
    }else{
      dirIdx=5;
    }
	} /*else{
		std::cout<< "something is weird in FeatureActiveParticlesTopological::checkMove "<< absX << " " << absY << " " << absZ << std::endl;
	} /* */

	return dirIdx;
}

#endif //LEMONADE_UPDATER_MOVES_MOVELOCALACTIVE_H
