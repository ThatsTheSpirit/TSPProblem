#pragma once
#include <vector>
#include <ctime>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <random>
#include <functional>
#include <fstream>
#include <chrono>


typedef std::vector<std::vector<int>> int_matrix;

class GamiltonPath
{
public:
	GamiltonPath();
	~GamiltonPath();

	void DoTest(int count);
	void ShowResults();
	void ShowBF();
	void ShowRS();
	void ShowG();
	void ShowSA();

	void SetTemperature(const int& temperature);
	void SetAlpha(const float& alpha);

	//int debugProperty ;//debug temp

private:
	int_matrix WeightGraph; //graph of all weights between each vertices //BF
	std::vector<int> Vertices; //vector of vertices //BF
	std::vector<int> VerticesRS; //vector of vertices //RS
	std::vector<int> VerticesG; //vector of vertices //Greedy
	std::vector<int> VerticesSA; //vector of vertices //SA

	std::vector<int> MinLengths; //BF
	std::vector<int> MinLengthsRS; //RS 
	std::vector<int> MinLengthsG; //G
	std::vector<int> MinLengthsSA; //SA

	float MiddleValueBF, MiddleValueRS, MiddleValueG, MiddleValueSA;
	float SigmaBF, SigmaRS, SigmaG, SigmaSA;
	float deltaQBF, deltaQRS, deltaQG, deltaQSA;
	float temperature;
	float alpha;

	std::pair<float, float> intervalBF;
	std::pair<float, float> intervalRS;
	std::pair<float, float> intervalG;
	std::pair<float, float> intervalSA;
		

	void writeToFile(const std::vector<int>& MinLengths, const std::string& name);

	void DoMinLenghtVectorBF();
	void DoMinLenghtVectorG();
	void DoMinLenghtVectorRS();
	void DoMinLengthVectorSA();

	int GetCost(const std::vector<int>& Vertices);

	void DoDeltaQ(float & deltaQ, const float & Sigma, const std::vector<int> & MinLengths);
	void DoRandVec();
	
	void DoMiddleValue(float & MiddleValue, const std::vector<int> & MinLengths);
	void DoSigma(float & Sigma, const float & MiddleValue, const std::vector<int> & MinLengths);
	void DoInterval(std::pair<float, float> & interval, const float & deltaQ, const float & MiddleValue);
};

