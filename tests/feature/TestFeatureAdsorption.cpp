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

#include <iostream>
#include <cstdio>

#include "gtest/gtest.h"

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureAdsorption.h>


class TestFeatureAdsorption: public ::testing::Test{
public:
  typedef LOKI_TYPELIST_1(FeatureAdsorption) Features;
  typedef ConfigureSystem<VectorInt3,Features,3> Config;
  typedef Ingredients<Config> IngredientsType;
  
  /*
  //redirect cout output
  virtual void SetUp(){
    originalBuffer=std::cout.rdbuf();
    std::cout.rdbuf(tempStream.rdbuf());
  };
  
  //restore original output
  virtual void TearDown(){
    std::cout.rdbuf(originalBuffer);
  };
  
private:
  std::streambuf* originalBuffer;
  std::ostringstream tempStream;
  */
};


TEST_F(TestFeatureAdsorption,Constructor)
{
  // instanz of ingredients should call constructor of FeatureAdsorption
  IngredientsType ingredients;
  
  // check the initial values of public accessibe variables
  EXPECT_DOUBLE_EQ(0.0,ingredients.getAdsorptionEnergyX());
  EXPECT_DOUBLE_EQ(0.0,ingredients.getAdsorptionEnergyY());
  EXPECT_DOUBLE_EQ(0.0,ingredients.getAdsorptionEnergyZ());
  
  EXPECT_FALSE(ingredients.getAdsorptionX());
  EXPECT_FALSE(ingredients.getAdsorptionY());
  EXPECT_FALSE(ingredients.getAdsorptionZ());
  
  // check if synchronize throws exception if no values are set
  EXPECT_ANY_THROW(ingredients.synchronize());

}

TEST_F(TestFeatureAdsorption,synchronize)
{
  IngredientsType ingredients;
  
  // if nothing is set, throw exception
  EXPECT_ANY_THROW(ingredients.synchronize());
  
  // set system parameters
  ingredients.setBoxX(16);
  ingredients.setBoxY(16);
  ingredients.setBoxZ(16);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  
  // throw still exception, no energy is set
  EXPECT_ANY_THROW(ingredients.synchronize());
  ingredients.setAdsorptionEnergyX(0.3);
  
  // throw still exception, box settings do not match walls
  EXPECT_ANY_THROW(ingredients.synchronize());
  ingredients.setPeriodicX(false);
  
  // now it should work
  EXPECT_NO_THROW(ingredients.synchronize());
  ingredients.setPeriodicY(false);
  ingredients.setAdsorptionEnergyY(0.4);
  ingredients.setPeriodicZ(false);
  ingredients.setAdsorptionEnergyZ(0.5);
  
  // now it should work again
  EXPECT_NO_THROW(ingredients.synchronize());
  ingredients.setPeriodicZ(true);
  
  // again, box does not match the adsorbing walls
  EXPECT_ANY_THROW(ingredients.synchronize());
  ingredients.setAdsorptionZ(false);
  
  // now it should work again
  EXPECT_NO_THROW(ingredients.synchronize());
  
}

   
TEST_F(TestFeatureAdsorption,ReadWriteRoutine)
{
  IngredientsType ingredients;
}
