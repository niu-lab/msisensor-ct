/*
 * homo.cpp for MSIsensor-ct 
 * Copyright (c) 2013 Beifang Niu && Kai Ye WUGSC All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include <iostream>
#include <sstream>
#include <bitset>
#include <omp.h>
#include <cmath>
#include <inttypes.h>
#include <fstream>

#include "homo.h"
#include "polyscan.h"
#include "somatic.h"
#include "xgboost/c_api.h"
#include "md5.hpp"

using namespace std;

extern Param param;
extern Param paramd;
extern PolyScan polyscan;
extern char homo_code[];
bool ctDNA; // 

HomoSite::HomoSite()
	: homoType(0)
	, length(0)
	, frontKmer(0)
	, endKmer(0)
	  , typeLen(0)
	  , chr("")
	  , transfer("")
	  , bases("")
	  , fbases("")
	, ebases("")
	, location(0)
	, lowcut(0)
	, highcut(0)
	, normalDis(NULL)
	, tumorDis(NULL)
	  , tumorDisNor(NULL)
	  ////////////
	, normalCov( 0 )
	, tumorCov( 0 )
	, withSufCov( false )
	, normalWithSufCov( false )
	, dif( 0.0 )
	, pValue( 0.0 )
	, comentropy( 0.0 )
	, somatic( false )
	  , withGenotype( false )   
{
	InitType();
};

HomoSite::~HomoSite() {
	// xxxxx
};

// transfer binary to string
void HomoSite::TransferString() {
	// assign memory
	//chr.resize(MAX_TRANSFER_LINE_LENGTH);
	std::stringstream ss;
	bit8_t tch = 0;
	char tchbuff[MAX_FLANK_REGION];
	for (int i=typeLen-1; i>=0; i--) {
		tch = 0;
		tch = homoType&3;
		tchbuff[i] = homo_code[tch];
		homoType = homoType>>2;
	}
	// load type
	for (int i=0; i<typeLen; i++) {
		ss<<tchbuff[i];
	}
	ss << '\t';
	for (int i=param.ContextLength-1; i>=0; i--) {
		tch = 0;
		tch = frontKmer&3;
		tchbuff[i] = homo_code[tch];
		frontKmer = frontKmer>>2;
	}
	// load front flank region
	for (int i=0; i<param.ContextLength; i++) {
		ss<<tchbuff[i];
	}
	ss << '\t';
	for (int i=param.ContextLength-1; i>=0; i--) {
		tch = 0;
		tch = endKmer&3;
		tchbuff[i] = homo_code[tch];
		endKmer = endKmer>>2;
	}
	// load end flank region
	for (int i=0; i<param.ContextLength; i++) {
		ss<<tchbuff[i];
	}
	transfer = ss.str();

	ss.clear();
	ss.str("");
}

// initial distribution
void HomoSite::InitialDis() {
	normalDis = new unsigned short *[polyscan.totalBamPairsNum];
	tumorDis  = new unsigned short *[polyscan.totalBamPairsNum];
	for (unsigned int j=0; j<polyscan.totalBamPairsNum; j++) {
		normalDis[j] = new unsigned short [paramd.s_dispots];
		tumorDis[j]  = new unsigned short [paramd.s_dispots];
		for (unsigned int k=0; k<paramd.s_dispots; k++) {
			normalDis[j][k] = 0;
			tumorDis[j][k]  = 0;
		}
	}
}

//initial tumor only distribution
void HomoSite::InitialTumorDis() {
	tumorDis  = new unsigned short *[polyscan.totalBamTumorsNum];
	for (unsigned int j=0; j<polyscan.totalBamTumorsNum; j++) {
		tumorDis[j]  = new unsigned short [paramd.s_dispots];
		for (unsigned int k=0; k<paramd.s_dispots; k++) {
			tumorDis[j][k]  = 0;
		}
	}
}

void HomoSite::InitialTumorDisNor() {
	tumorDisNor = new float [paramd.t_dispots];
	for (unsigned int n=0; n<paramd.t_dispots; n++) {
		tumorDisNor[n] = 0.0;
	}
}

// Out distribution
void HomoSite::OutputDis() {
	for (unsigned int j=0; j<polyscan.totalBamPairsNum; j++) {
		for (unsigned int k=0; k<paramd.s_dispots; k++) {
			std::cout<<normalDis[j][k];
		}
		for (unsigned int k=0; k<paramd.s_dispots; k++) {
			std::cout<<tumorDis[j][k];
		}
	}
}

//Out tumor only distribution
void HomoSite::OutputTumorDis() {
	for (unsigned int j=0; j<polyscan.totalBamTumorsNum; j++) {
		for (unsigned int k=0; k<paramd.s_dispots; k++) {
			std::cout<<tumorDis[j][k];
		}
	}
}

// Pourout distribution
void HomoSite::PouroutDis(Sample &sample) {
	sample.outputDistribution << chr << " "
		<< location << " " 
		<< fbases << " "
		<< length << "["
		<< bases  << "] "
		<< ebases <<"\n";

	for (unsigned int j=0; j<polyscan.totalBamPairsNum; j++) {
		sample.outputDistribution << "N: ";
		for (unsigned int k=0; k<paramd.s_dispots; k++) {
			sample.outputDistribution << normalDis[j][k] << " ";
		}
		sample.outputDistribution << "\nT: ";
		for (unsigned int k=0; k<paramd.s_dispots; k++) {
			sample.outputDistribution << tumorDis[j][k] << " ";
		}

	}
	sample.outputDistribution << "\n"; 
}

//Pour tumor only distribution
void HomoSite::PourTumoroutDis(Sample &sample) {
	sample.outputDistribution << chr << " "
		<< location << " "
		<< fbases << " "
		<< length << "["
		<< bases  << "] "
		<< ebases <<"\n";

	for (unsigned int j=0; j<polyscan.totalBamTumorsNum; j++) {
		sample.outputDistribution << "T: ";
		for (unsigned int k=0; k<paramd.s_dispots; k++) {
			sample.outputDistribution << tumorDis[j][k] << " ";
		}   
	}
	sample.outputDistribution << "\n";
}

// initial bools
void HomoSite::BoolsInitial() {
	withSufCov = normalWithSufCov = somatic = withGenotype = false;
}

// genotyping analyis
void HomoSite::DisGenotyping(Sample &sample) {
	/////// genotyping //////////
	// BoolsInitial();
	// Total sites
	////////////////////////////
	sample.numberOftotalSites++;
	///////////////////////////
	bool reportSomatic, reportGermline;
	reportSomatic = reportGermline = false;
	// update all 
	normalCov = tumorCov = 0;
	for (unsigned int j=0; j<polyscan.totalBamPairsNum; j++) {
		for (unsigned int k=0; k<paramd.s_dispots; k++) {
			normalCov += normalDis[j][k];
			tumorCov  += tumorDis[j][k];
		}
	}

	normalWithSufCov = (normalCov >= paramd.covCutoff) ? true : false;

	if ( normalWithSufCov && (tumorCov >= paramd.covCutoff) ) {
		withSufCov = true;
		dif = DistanceBetweenTwo( normalDis[0], tumorDis[0] );
		pValue = X2BetweenTwo( normalDis[0], tumorDis[0], paramd.s_dispots );
		// if (pValue < 0.001) somatic = true;
		// add one for FDR
		somatic = true;

	} else {

		withSufCov = false;
		dif = -1.0;
		pValue = 1;
	}
	// compute genotype
	ComputeGenotype( normalDis[0] );
	///////////////////////////////////
	unsigned PairIndex = 0;
	if ( withSufCov  ) {
		sample.numberOfDataPoints++;
		//sample.outputPSomatic << log10( pValue ) << std::endl;
		//if ( somatic ) sample.numberOfMsiDataPoints++;
	}

	if ( withGenotype ) reportGermline = true;
	if ( somatic ) reportSomatic = true;

	if ( reportSomatic ) {
		/*
		   sample.outputSomatic << chr << "\t" 
		   << location << "\t" 
		   << fbases << "\t" 
		   << length << "\t"
		   << bases << "\t" 
		   << ebases;
		   sample.outputSomatic << "\t" << std::fixed << dif << "\t"<< std::fixed << pValue;
		   sample.outputSomatic << std::endl;
		 */
		// record for FDR analysis 
		// instead of report it directly
		SomaticSite onessite;
		onessite.chr = chr;
		onessite.location = location;
		onessite.length = length;
		onessite.fbases = fbases;
		onessite.ebases = ebases;
		onessite.bases = bases;
		onessite.diff = dif;
		onessite.pValue = pValue;

		sample.totalSomaticSites.push_back( onessite );
	}

	if ( reportGermline ) {
		sample.outputGermline << chr << "\t" 
			<< location << "\t" 
			<< fbases << "\t" 
			<< length << "\t" 
			<< bases << "\t" 
			<< ebases;
		sample.outputGermline << "\t" << genotype[0] << "|" << genotype[1];
		sample.outputGermline << std::endl;
	}

}

// tumor somatic analyis
void HomoSite::DisTumorSomatic(Sample &sample) {
	sample.numberOftotalSites++;
	bool reportSomatic;
	reportSomatic = false;
	tumorCov = 0;
	for (unsigned int j=0; j<polyscan.totalBamTumorsNum; j++) {
		for (unsigned int k=0; k<paramd.s_dispots; k++) {
			tumorCov  += tumorDis[j][k];
		}
	}

	if (tumorCov >= paramd.covCutoff) {
		withSufCov = true;
		comentropy = Comentropy( tumorDis[0], paramd.s_dispots);
	} else {
		withSufCov = false;
		comentropy = 0;
	}
	//if (comentropy >= paramd.comentropyThreshold) {
	if (withSufCov and comentropy >= paramd.comentropyThreshold) { 
		reportSomatic = true;
		sample.numberOfMsiDataPoints ++;
	}
	if (withSufCov) sample.numberOfDataPoints++;
	if (reportSomatic) {
		sample.outputSomatic << chr << "\t"
			<< location << "\t"
			<< fbases << "\t"
			<< length << "\t"
			<< bases << "\t"
			<< ebases;
		sample.outputSomatic << "\t" << std::fixed << comentropy;
		sample.outputSomatic << std::endl;
		SomaticSite onessite;
		onessite.chr = chr;
		onessite.location = location;
		onessite.length = length;
		onessite.fbases = fbases;
		onessite.ebases = ebases;
		onessite.bases = bases;

		sample.totalSomaticSites.push_back( onessite );
	}
}

std::string HomoSite::DecIntToHexStr(long long num) {
	std::string str;
	long long Temp = num / 16;
	int left = num % 16;
	if (Temp > 0)
		str += DecIntToHexStr(Temp);
	if (left < 10)
		str += (left + '0');
	else
		str += ('a' + left - 10);
	return str;
}

std::string HomoSite::encode(const std::string &fileName) {
	unsigned char decrypt[16];
	std::string decryFileName;
	unsigned char* encrypt = (unsigned char*)fileName.c_str();
	MD5_CTX md5;
	MD5Init(&md5);
	MD5Update(&md5, encrypt, (int)strlen((char *)encrypt));
	MD5Final(&md5, decrypt);
	for (int i = 0; i<16; i++){
		string str = DecIntToHexStr((long long)decrypt[i]);
		while(str.size() < 2)
			str = "0" + str;
		decryFileName += str;
	}
	return decryFileName;
}

void HomoSite::DisTumorFeature(Sample &sample, const std::string &modelPre) {
	fstream _file;
	bool reportSomatic;
	std::string siteLoc;
	std::string ModelFile;
	std::string entryFileName;
	const float *m_TestResults2;
	reportSomatic = false;
	tumorCov = 0;
	sample.numberOftotalSites++;
	for (unsigned int j=0; j<polyscan.totalBamTumorsNum; j++) {
		for (unsigned int k=0; k<paramd.s_dispots; k++) {
			tumorCov  += tumorDis[j][k];

		}
	}
	if (tumorCov >= paramd.covCutoff){
		withSufCov = true;
		Normalization( tumorDis[0], paramd.t_dispots,length);
		stringstream tmp;
		tmp<<location;
		string tmp_string;
		tmp>>tmp_string;
		siteLoc = chr + "_" + tmp_string;
		entryFileName = encode(siteLoc);
		ModelFile =  modelPre + "/" + entryFileName;
		_file.open(ModelFile.c_str());
		if (_file) {
			sample.numberOfDataPoints++;
			BoosterHandle handle;
			XGBoosterCreate(NULL, 0, &handle);
			int x = XGBoosterLoadModel(handle, ModelFile.c_str());
			if (x != 0) {
				cout << "Failed Loading Model" << endl;
				exit(-1);
			}
			DMatrixHandle dHandle;
			//for(int i=0; i<paramd.s_dispots; i++){
             			//cout<< tumorDisNor[i] << " " ;
           			// }			 
            			//cout<< "/n" << endl;
            x = XGDMatrixCreateFromMat(tumorDisNor, 1, paramd.t_dispots-1, 0, &dHandle);
			if (x != 0)
			{
				cout << "Failed Creating DMatrix" << endl;
				exit(-1);
			}
			bst_ulong out2;
			x = XGBoosterPredict(handle, dHandle, 0, 0, &out2, &m_TestResults2);
			if (x != 0) {
				cout << "Failed Predicting" << endl;
				exit(-1);
			}
			if(m_TestResults2[0] > 0.5) {
				reportSomatic = true;
				sample.numberOfMsiDataPoints ++;
			}
			XGBoosterFree(handle);
		}
		_file.close();
		if (reportSomatic) {
			sample.outputSomatic << chr << "\t"
				<< location << "\t"
				<< fbases << "\t"
				<< length << "\t"
				<< bases << "\t"
				<< ebases;
			sample.outputSomatic << "\t" << std::fixed << m_TestResults2[0];
			sample.outputSomatic << std::endl;
			SomaticSite onessite;
			onessite.chr = chr;
			onessite.location = location;
			onessite.length = length;
			onessite.fbases = fbases;
			onessite.ebases = ebases;
			onessite.bases = bases;

			sample.totalSomaticSites.push_back( onessite );
		}
	}
}

// distance
double HomoSite::DistanceBetweenTwo( unsigned short * FirstOriginal, unsigned short * SecondOriginal ) {
	double SmallDouble = 0.0000000001;
	// declare 
	double *Min, *Max, *FirstNormalized, *SecondNormalized, AreaMin, AreaMax;
	unsigned int dispots, sumFirst, sumSecond;
	dispots = paramd.s_dispots;
	FirstNormalized   = new double [dispots + 2];
	SecondNormalized  = new double [dispots + 2];
	sumFirst = sumSecond = 0;
	// sum 
	for (unsigned i = 0; i< dispots; i++) {
		sumFirst  += FirstOriginal[i];
		sumSecond += SecondOriginal[i];
	}
	// normalization
	for (unsigned i = 0; i < dispots; i++) {
		if (FirstOriginal[i] < 1) FirstNormalized[i]  = 0.0;
		else FirstNormalized[i]  = (FirstOriginal[i] / (double)sumFirst);
		if (SecondOriginal[i] < 1) SecondNormalized[i] = 0.0;
		else SecondNormalized[i] = SecondOriginal[i] / (double)sumSecond;
	}
	FirstNormalized[dispots]  = SecondNormalized[dispots] = 0.0;
	Min  = new double [dispots + 2];
	Max  = new double [dispots + 2];
	// initialization 
	for (unsigned i=0; i<=dispots; i++) Min[i] = Max[i] = 0;
	// get min and max per position
	for (unsigned i=0; i<=dispots; i++) {
		if (FirstNormalized[i] <= SecondNormalized[i]) {
			Min[i] = FirstNormalized[i];
			Max[i] = SecondNormalized[i];
		} else {
			Min[i] = SecondNormalized[i];
			Max[i] = FirstNormalized[i];
		}
		if (Min[i] < SmallDouble) Min[i] = 0.0;
		if (Max[i] < SmallDouble) Max[i] = 0.0;
	}

	// caculate area
	AreaMin =  AreaMax = 0.0;
	for (unsigned i=0; i< dispots; i++) {
		AreaMin += (Min[i] + Min[i + 1]);
		AreaMax += (Max[i] + Max[i + 1]);
	}
	if (AreaMax < AreaMin) {
		std::cerr << "something is wrong with distance calculation " << AreaMax << " " << AreaMin << std::endl;
		for (unsigned index = 0; index < 100; index++) {
			std::cerr << index << "\t" 
				<< FirstOriginal[index] << "\t" 
				<< SecondOriginal[index] << "\t" 
				<< sumFirst << "\t" 
				<< sumSecond << "\t" 
				<< FirstNormalized[index] << "\t" 
				<< SecondNormalized[index] << "\t" 
				<< Min[index] << "\t" 
				<< Max[index] 
				<< std::endl;
		}
	}
	delete [] Min;
	delete [] Max;

	delete [] FirstNormalized;
	delete [] SecondNormalized;

	return (AreaMax - AreaMin)/AreaMax;
}


//comentropy
double HomoSite::Comentropy( unsigned short * tumorDis, unsigned int dispots ) {
	double sum = 0;
	double comentropy = 0.0;
	double number = 0.0;
	for (int i = 0; i < dispots; i++){
		if (tumorDis[i] <3){
			tumorDis[i] = 0;
		}
		sum += tumorDis[i];
	}
	if ( sum != 0 ) {
		for( int j = 0; j < dispots; j++){
			if( tumorDis[j] != 0 ){
				number += 1;
				comentropy -= (tumorDis[j]/sum)*log(tumorDis[j]/sum);
			}
		}
	}
	if ( number != 0 ) comentropy = comentropy/number;
	return comentropy;
}


void HomoSite::Normalization( unsigned short * tumorDis, unsigned int dispots , unsigned int length ) {
	if (ctDNA)
	{
		//ofstream _debug("debug.2020"); // add for debug
		unsigned short * local_tumorDis = new unsigned short[dispots] ;
		for (int i = 0; i < dispots; i++){
                        local_tumorDis[i] = tumorDis[i];
		}
        
          
              	local_tumorDis[length-1] = 0;
                cout << "quick" << " " ;
                cout << length << " ";
		cout << "/n" << endl;
		//cout<< "quick" << " ";
		for(int i=0;i<dispots;i++){
		cout << local_tumorDis[i] << " ";
		}
		cout << "/n" << endl;
        
		double sum = 0;
		for (int i = 0; i < dispots; i++){
			sum += local_tumorDis[i];
		}
		for (int j = 0; j < dispots; j++){
			if(abs(sum) <= 1e-15){
			tumorDisNor[j] = local_tumorDis[j];
			}
			else
			{tumorDisNor[j] = local_tumorDis[j]/sum;}
		}
		delete []local_tumorDis;

	}
	else
	{
		double sum = 0;
		for (int i = 0; i < dispots; i++){
			sum += tumorDis[i];
		}
		for (int j = 0; j < dispots; j++){
			tumorDisNor[j] = tumorDis[j]/sum;
		}
	}
}

void HomoSite::ComputeGenotype( unsigned short * NormalReadCount ) {
	unsigned int Offset, CoverageCutoff, first, second, Sum;
	Offset = 1; CoverageCutoff = paramd.covCutoff;
	first = second = Sum = 0;

	// find the largest number 
	for (unsigned int pos_index = 0; pos_index < paramd.s_dispots; pos_index++) { 
		// NormalReadCount
		Sum += NormalReadCount[pos_index];
		if (NormalReadCount[pos_index] > NormalReadCount[first]) first = pos_index;
	}

	if (Sum < CoverageCutoff) {
		genotype[0] = genotype[1] = -1;
		withGenotype = false;
		return;
	}

	if (first == 0) second = 1;
	for (unsigned int pos_index = 0; pos_index < paramd.s_dispots; pos_index++) {
		if (pos_index == first) continue; 
		if (NormalReadCount[pos_index] > NormalReadCount[second]) second = pos_index;
	}
	float first_ratio  = NormalReadCount[first]  / (float) Sum;
	float second_ratio = NormalReadCount[second] / (float) Sum;

	if (first_ratio > 0.7) {
		genotype[0] = genotype[1]  = first + Offset;
		withGenotype = true; 
	}
	else if (first_ratio > 0.3 && second_ratio < first_ratio * 0.66) {
		genotype[0] =  genotype[1]  = first + Offset;
		withGenotype = true; 
	}
	else if (first_ratio > 0.3 && second_ratio >= first_ratio * 0.66) {
		genotype[0] = first + Offset; 
		genotype[1] = second + Offset; 
		withGenotype = true; 
	}

	return;
} 

// release memory
void HomoSite::ReleaseMemory() {
	for (unsigned int k=0; k<polyscan.totalBamPairsNum; k++) {
		delete [] normalDis[k];
		delete []  tumorDis[k]; 
	}
	delete [] normalDis;
	delete [] tumorDis;
}

//release tumor memory 
void HomoSite::ReleaseTumorMemory() {
	for (unsigned int k=0; k<polyscan.totalBamTumorsNum; k++) {
		delete []  tumorDis[k]; 
	}
	delete [] tumorDis;
}

void HomoSite::ReleaseTumorMemoryNor() {
	for (unsigned int k=0; k<polyscan.totalBamTumorsNum; k++) {
		delete []  tumorDis[k]; 
	}
	delete [] tumorDis;
	delete [] tumorDisNor;
}
