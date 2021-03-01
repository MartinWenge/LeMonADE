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

/*****************************************************************************/
/**
 * @file
 * @brief test for Feature Tanglotron
 * @author Martin Wengenmayr
 * @date 01.03.2021
 * */
/*****************************************************************************/
#include "gtest/gtest.h"


#include <iostream>

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>

#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveLocalScDiag.h>

#include <LeMonADE/feature/FeatureTanglotron.h>


class TestFeatureTanglotron: public ::testing::Test{
public:
  typedef LOKI_TYPELIST_2(FeatureMoleculesIO,FeatureTanglotron) Features;
  typedef ConfigureSystem<VectorInt3,Features> ConfigType;
  typedef Ingredients < ConfigType> IngredientsType;
  IngredientsType ingredients;

  //dummy move class used to check response to unknown move type
  class UnknownMove:public MoveBase
  {
	  public:
		template<class IngredientsType> bool check(IngredientsType& ingredients) const
		{
		  return ingredients.checkMove(ingredients,*this);
		}

		template<class IngredientsType> void apply(IngredientsType& ingredients)
		{
		  ingredients.applyMove(ingredients,*this);
		}

		template <class IngredientsType> void init(const IngredientsType& ingredients){};
  };
  
  /* suppress cout output for better readability -->un-/comment here: */
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
  /* ** */
};


TEST_F(TestFeatureTanglotron,Basics){
  // test default constructor
  FeatureTanglotron feature;
  EXPECT_TRUE(feature.getCalculatePotential());
  EXPECT_DOUBLE_EQ(0.0,feature.getTorque());
  EXPECT_TRUE(feature.getTanglotronMotorMap().empty());
  // modify system variables
  feature.setTorque(3.14);
  feature.setCalculatePotential(false);
  EXPECT_FALSE(feature.getCalculatePotential());
  EXPECT_DOUBLE_EQ(3.14,feature.getTorque());
  
  //test standard constructor in ingredients
  EXPECT_TRUE(ingredients.getCalculatePotential());
  EXPECT_DOUBLE_EQ(0.0,ingredients.getTorque());

}

TEST_F(TestFeatureTanglotron, SystemInformation){
  //setup a small system
  ingredients.setBoxX(16);
  ingredients.setBoxY(16);
  ingredients.setBoxZ(16);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.setTorque(3.14);
  ingredients.modifyBondset().addBFMclassicBondset();

  ingredients.modifyMolecules().addMonomer(0,0,8);
  // check default values
  EXPECT_EQ(-1, ingredients.getMolecules()[0].getTanglotronID());
  EXPECT_EQ(-1, ingredients.getMolecules()[0].getTanglotronType());
  EXPECT_FALSE(ingredients.getMolecules()[0].getIsTanglotron());
  ingredients.modifyMolecules()[0].setIsTanglotron(true);
  ingredients.modifyMolecules()[0].setTanglotronType(TanglotronAttributeTag::tanglotronType::stator);
  ingredients.modifyMolecules()[0].setTanglotronID(0);
  ingredients.modifyMolecules().addMonomer(2,1,8);
  ingredients.modifyMolecules()[1].setIsTanglotron(true);
  ingredients.modifyMolecules()[1].setTanglotronType(TanglotronAttributeTag::tanglotronType::rotorA);
  ingredients.modifyMolecules()[1].setTanglotronID(0);
  ingredients.modifyMolecules().addMonomer(2,-1,8);
  ingredients.modifyMolecules()[2].setIsTanglotron(true);
  ingredients.modifyMolecules()[2].setTanglotronType(TanglotronAttributeTag::tanglotronType::rotorB);
  ingredients.modifyMolecules()[2].setTanglotronID(0);
  ingredients.modifyMolecules().connect(0,1);
  ingredients.modifyMolecules().connect(0,2);
  ingredients.modifyMolecules().connect(1,2);
  ingredients.modifyMolecules().addMonomer(-2,1,8);
  ingredients.modifyMolecules()[3].setIsTanglotron(true);
  ingredients.modifyMolecules()[3].setTanglotronType(TanglotronAttributeTag::tanglotronType::rotorC);
  ingredients.modifyMolecules()[3].setTanglotronID(0);
  ingredients.modifyMolecules().addMonomer(-2,-1,8);
  ingredients.modifyMolecules()[4].setIsTanglotron(true);
  ingredients.modifyMolecules()[4].setTanglotronType(TanglotronAttributeTag::tanglotronType::rotorD);
  ingredients.modifyMolecules()[4].setTanglotronID(0);
  ingredients.modifyMolecules().connect(0,3);
  ingredients.modifyMolecules().connect(0,4);
  ingredients.modifyMolecules().connect(3,4);

  // check the map before sync
  EXPECT_TRUE(ingredients.getTanglotronMotorMap().empty());
  // snchronize
  ingredients.synchronize();
  // check the map after sync
  EXPECT_EQ(1, ingredients.getTanglotronMotorMap().size());
  EXPECT_EQ(5, ingredients.getTanglotronMotorMap().at(0).size());
  EXPECT_EQ(0, ingredients.getTanglotronMotorMap().at(0).at(0));
  EXPECT_EQ(1, ingredients.getTanglotronMotorMap().at(0).at(1));
  EXPECT_EQ(2, ingredients.getTanglotronMotorMap().at(0).at(2));
  EXPECT_EQ(3, ingredients.getTanglotronMotorMap().at(0).at(3));
  EXPECT_EQ(4, ingredients.getTanglotronMotorMap().at(0).at(4));

  // check base move
  // UnknownMove basemove;
  // basemove.init(ingredients);
  // EXPECT_TRUE(basemove.check(ingredients));
  // //should change nothing
  // basemove.apply(ingredients);
  // EXPECT_EQ(1.0, basemove.getProbability());

  // check move local sc
  MoveLocalSc localmove;
  // possible stator move in z direction
  localmove.init(ingredients, 0, VectorInt3(0,0,1));
  EXPECT_TRUE(localmove.check(ingredients));
  EXPECT_EQ(1.0, localmove.getProbability());

  // move hindered by bond vector constraints
  localmove.init(ingredients, 0, VectorInt3(1,0,0));
  EXPECT_FALSE(localmove.check(ingredients));
  EXPECT_EQ(1.0, localmove.getProbability());

  localmove.init(ingredients, 1, VectorInt3(0,0,1));
  EXPECT_TRUE(localmove.check(ingredients));
  EXPECT_NE(1.0, localmove.getProbability());
}

TEST_F(TestFeatureTanglotron, ReadWrite){
  //setup a small system
  ingredients.setBoxX(16);
  ingredients.setBoxY(16);
  ingredients.setBoxZ(32);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.setTorque(3.14);
  ingredients.modifyBondset().addBFMclassicBondset();

  ingredients.modifyMolecules().addMonomer(0,0,8);
  ingredients.modifyMolecules()[0].setIsTanglotron(true);
  ingredients.modifyMolecules()[0].setTanglotronType(TanglotronAttributeTag::tanglotronType::stator);
  ingredients.modifyMolecules()[0].setTanglotronID(0);
  ingredients.modifyMolecules().addMonomer(2,1,8);
  ingredients.modifyMolecules()[1].setIsTanglotron(true);
  ingredients.modifyMolecules()[1].setTanglotronType(TanglotronAttributeTag::tanglotronType::rotorA);
  ingredients.modifyMolecules()[1].setTanglotronID(0);
  ingredients.modifyMolecules().addMonomer(2,-1,8);
  ingredients.modifyMolecules()[2].setIsTanglotron(true);
  ingredients.modifyMolecules()[2].setTanglotronType(TanglotronAttributeTag::tanglotronType::rotorB);
  ingredients.modifyMolecules()[2].setTanglotronID(0);
  ingredients.modifyMolecules().connect(0,1);
  ingredients.modifyMolecules().connect(0,2);
  ingredients.modifyMolecules().connect(1,2);
  ingredients.modifyMolecules().addMonomer(-2,1,8);
  ingredients.modifyMolecules()[3].setIsTanglotron(true);
  ingredients.modifyMolecules()[3].setTanglotronType(TanglotronAttributeTag::tanglotronType::rotorC);
  ingredients.modifyMolecules()[3].setTanglotronID(0);
  ingredients.modifyMolecules().addMonomer(-2,-1,8);
  ingredients.modifyMolecules()[4].setIsTanglotron(true);
  ingredients.modifyMolecules()[4].setTanglotronType(TanglotronAttributeTag::tanglotronType::rotorD);
  ingredients.modifyMolecules()[4].setTanglotronID(0);
  ingredients.modifyMolecules().connect(0,3);
  ingredients.modifyMolecules().connect(0,4);
  ingredients.modifyMolecules().connect(3,4);

  ingredients.modifyMolecules().addMonomer(0,0,0);

  // snchronize
  ingredients.synchronize();

  //write to file and read back in
  AnalyzerWriteBfmFile<IngredientsType> anaWrite("test_feature_tanglotron.bfm", ingredients, AnalyzerWriteBfmFile<IngredientsType>::OVERWRITE);
  anaWrite.initialize();
  anaWrite.execute();
  anaWrite.cleanup();

  IngredientsType readInIngredients;
  UpdaterReadBfmFile<IngredientsType> updRead("test_feature_tanglotron.bfm",readInIngredients,UpdaterReadBfmFile<IngredientsType>::READ_LAST_CONFIG_SAVE);
  updRead.initialize();
  updRead.execute();
  updRead.cleanup();

  remove("test_feature_tanglotron.bfm");

  EXPECT_DOUBLE_EQ(3.14, readInIngredients.getTorque());
  EXPECT_EQ(6, readInIngredients.getMolecules().size());
  EXPECT_EQ(16, readInIngredients.getBoxX());
  EXPECT_EQ(16, readInIngredients.getBoxY());
  EXPECT_EQ(32, readInIngredients.getBoxZ());
  EXPECT_TRUE(readInIngredients.isPeriodicX());
  EXPECT_TRUE(readInIngredients.isPeriodicY());
  EXPECT_TRUE(readInIngredients.isPeriodicZ());

  EXPECT_TRUE(readInIngredients.getMolecules()[0].getIsTanglotron());
  EXPECT_EQ(0,readInIngredients.getMolecules()[0].getTanglotronID());
  EXPECT_EQ(0,readInIngredients.getMolecules()[0].getTanglotronType());
  EXPECT_TRUE(readInIngredients.getMolecules()[1].getIsTanglotron());
  EXPECT_EQ(0,readInIngredients.getMolecules()[1].getTanglotronID());
  EXPECT_EQ(1,readInIngredients.getMolecules()[1].getTanglotronType());
  EXPECT_TRUE(readInIngredients.getMolecules()[2].getIsTanglotron());
  EXPECT_EQ(0,readInIngredients.getMolecules()[2].getTanglotronID());
  EXPECT_EQ(2,readInIngredients.getMolecules()[2].getTanglotronType());
  EXPECT_TRUE(readInIngredients.getMolecules()[3].getIsTanglotron());
  EXPECT_EQ(0,readInIngredients.getMolecules()[3].getTanglotronID());
  EXPECT_EQ(3,readInIngredients.getMolecules()[3].getTanglotronType());
  EXPECT_TRUE(readInIngredients.getMolecules()[4].getIsTanglotron());
  EXPECT_EQ(0,readInIngredients.getMolecules()[4].getTanglotronID());
  EXPECT_EQ(4,readInIngredients.getMolecules()[4].getTanglotronType());

  EXPECT_FALSE(readInIngredients.getMolecules()[5].getIsTanglotron());
  EXPECT_EQ(-1,readInIngredients.getMolecules()[5].getTanglotronID());
  EXPECT_EQ(-1,readInIngredients.getMolecules()[5].getTanglotronType());

  EXPECT_EQ(VectorInt3(0,0,8),readInIngredients.getMolecules()[0].getVector3D());
  EXPECT_EQ(VectorInt3(2,1,8),readInIngredients.getMolecules()[1].getVector3D());
  EXPECT_EQ(VectorInt3(2,-1,8),readInIngredients.getMolecules()[2].getVector3D());
  EXPECT_EQ(VectorInt3(-2,1,8),readInIngredients.getMolecules()[3].getVector3D());
  EXPECT_EQ(VectorInt3(-2,-1,8),readInIngredients.getMolecules()[4].getVector3D());
  EXPECT_EQ(VectorInt3(0,0,0),readInIngredients.getMolecules()[5].getVector3D());

  EXPECT_TRUE(readInIngredients.getMolecules().areConnected(0,1));
  EXPECT_TRUE(readInIngredients.getMolecules().areConnected(0,2));
  EXPECT_TRUE(readInIngredients.getMolecules().areConnected(0,3));
  EXPECT_TRUE(readInIngredients.getMolecules().areConnected(0,4));
  EXPECT_TRUE(readInIngredients.getMolecules().areConnected(1,2));
  EXPECT_TRUE(readInIngredients.getMolecules().areConnected(3,4));

}
// TEST_F(TestFeatureTanglotron,fileReadWrite){
//   IngredientsType ingredientsWrite;
//   IngredientsType ingredientsRead;
// 
//   ingredientsWrite.setBoxX(64);
//   ingredientsWrite.setPeriodicX(true);
//   ingredientsWrite.setBoxY(64);
//   ingredientsWrite.setPeriodicY(true);
//   ingredientsWrite.setBoxZ(64);
//   ingredientsWrite.setPeriodicZ(true);
//   ingredientsWrite.setSpringConstant(1.2);
//   ingredientsWrite.setEquilibriumLength(4);
//   ingredientsWrite.modifyMolecules().addMonomer(0,0,0);
//   ingredientsWrite.modifyMolecules()[0].setMonomerGroupTag(1);
//   ingredientsWrite.modifyMolecules().addMonomer(2,0,0);
//   ingredientsWrite.modifyMolecules()[1].setMonomerGroupTag(1);
//   ingredientsWrite.modifyMolecules().addMonomer(4,0,0);
//   ingredientsWrite.modifyMolecules()[2].setMonomerGroupTag(2);
//   ingredientsWrite.modifyMolecules().addMonomer(6,0,0);
//   ingredientsWrite.modifyMolecules()[3].setMonomerGroupTag(2);
//   ingredientsWrite.modifyMolecules().addMonomer(8,0,0);
// 
//   EXPECT_EQ(1,ingredientsWrite.getMolecules()[0].getMonomerGroupTag());
//   EXPECT_EQ(1,ingredientsWrite.getMolecules()[1].getMonomerGroupTag());
//   EXPECT_EQ(2,ingredientsWrite.getMolecules()[2].getMonomerGroupTag());
//   EXPECT_EQ(2,ingredientsWrite.getMolecules()[3].getMonomerGroupTag());
//   EXPECT_EQ(0,ingredientsWrite.getMolecules()[4].getMonomerGroupTag());
// 
//   //write to file and read back in
//   AnalyzerWriteBfmFile<IngredientsType> anaWrite("tests/springPotentialTestOut.test",ingredientsWrite,AnalyzerWriteBfmFile<IngredientsType>::NEWFILE);
//   anaWrite.initialize();
//   anaWrite.execute();
// 
//   FileImport<IngredientsType> infile ("tests/springPotentialTestOut.test",ingredientsRead);
// 
//   //scan file for !mcs and read-in first frame
//   infile.initialize();
// 
//   EXPECT_EQ(ingredientsRead.getMolecules()[0].getMonomerGroupTag(),1);
//   EXPECT_EQ(ingredientsRead.getMolecules()[1].getMonomerGroupTag(),1);
//   EXPECT_EQ(ingredientsRead.getMolecules()[2].getMonomerGroupTag(),2);
//   EXPECT_EQ(ingredientsRead.getMolecules()[3].getMonomerGroupTag(),2);
//   EXPECT_EQ(ingredientsRead.getMolecules()[4].getMonomerGroupTag(),0); /*has default value*/
// 
//   EXPECT_DOUBLE_EQ(1.2,ingredientsRead.getSpringConstant());
//   EXPECT_DOUBLE_EQ(4.0,ingredientsRead.getEquilibriumLength());
// 
//     //remove the temporary file
//   EXPECT_EQ(0,remove("tests/springPotentialTestOut.test"));
// }
