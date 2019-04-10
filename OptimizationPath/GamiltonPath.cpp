#include "GamiltonPath.h"



GamiltonPath::GamiltonPath()
{
	SigmaBF = SigmaRS = SigmaG = SigmaSA = MiddleValueSA = MiddleValueBF =
		MiddleValueRS = MiddleValueG = debugProperty = deltaQBF = deltaQRS =
		deltaQG = deltaQSA = 0;

	temperature = 1000.0;
	alpha = 0.95;
		
}


void GamiltonPath::DoTest(int count)
{
	while (count--)
	{
		DoRandVec();
		DoMinLenghtVectorBF();//1
		DoMinLenghtVectorRS();//2
		DoMinLenghtVectorG();//3
		DoMinLengthVectorSA();//4
		
	}

	DoMiddleValue(MiddleValueBF, MinLengths);
	DoMiddleValue(MiddleValueRS, MinLengthsRS);
	DoMiddleValue(MiddleValueG, MinLengthsG);
	DoMiddleValue(MiddleValueSA, MinLengthsSA);

	DoSigma(SigmaBF, MiddleValueBF, MinLengths);
	DoSigma(SigmaRS, MiddleValueRS, MinLengthsRS);
	DoSigma(SigmaG, MiddleValueG, MinLengthsG);
	DoSigma(SigmaSA, MiddleValueSA, MinLengthsSA);

	DoDeltaQ(deltaQBF, SigmaBF, MinLengths);
	DoDeltaQ(deltaQRS, SigmaRS, MinLengthsRS);
	DoDeltaQ(deltaQG, SigmaG, MinLengthsG);
	DoDeltaQ(deltaQSA, SigmaSA, MinLengthsSA);

	DoInterval(intervalBF,deltaQBF,MiddleValueBF);
	DoInterval(intervalRS, deltaQRS, MiddleValueRS);
	DoInterval(intervalG, deltaQG, MiddleValueG);
	DoInterval(intervalSA, deltaQSA, MiddleValueSA);

	//writeToFile(MinLengths, "BF1000.txt");
	//writeToFile(MinLengthsRS, "RS1000.txt");
	//writeToFile(MinLengthsG, "G1000.txt");
	//writeToFile(MinLengthsSA, "SA1000.txt");
}

void GamiltonPath::ShowBF()
{
	std::cout << "\t\t\t\tBF Method:\n\n";
	std::cout << "MiddleLength = " << MiddleValueBF << "\tSigma = " << SigmaBF
		<< "\tDeltaQ = " << deltaQBF << std::endl
		<< "Confedience interval = {" << intervalBF.first << ", " << intervalBF.second << "}";
}

void GamiltonPath::ShowRS()
{
	std::cout << std::endl << std::endl << "\t\t\t\tRS Method:\n\n";
	std::cout << "MiddleLength = " << MiddleValueRS << "\tSigma = " << SigmaRS
		<< "\tDeltaQ = " << deltaQRS << std::endl
		<< "Confedience interval = {" << intervalRS.first << ", " << intervalRS.second << "}";
}

void GamiltonPath::ShowG()
{
	std::cout << std::endl << std::endl << "\t\t\t\tG Method:\n\n";
	std::cout << "MiddleLength = " << MiddleValueG << "\tSigma = " << SigmaG
		<< "\tDeltaQ = " << deltaQG << std::endl
		<< "Confedience interval = {" << intervalG.first << ", " << intervalG.second << "}";
}

void GamiltonPath::ShowSA()
{
	std::cout << std::endl << std::endl << "\t\t\t\tSA Method:\n\n";
	std::cout << "MiddleLength = " << MiddleValueSA << "\tSigma = " << SigmaSA
		<< "\tDeltaQ = " << deltaQSA << std::endl
		<< "Confedience interval = {" << intervalSA.first << ", " << intervalSA.second << "}";
}


void GamiltonPath::SetTemperature(const int & temperature)
{
	this->temperature = temperature;
}

void GamiltonPath::SetAlpha(const float & alpha)
{
	this->alpha = alpha;
}

void GamiltonPath::ShowResults()
{
	ShowBF();
	ShowRS();
	ShowG();
	ShowSA();
	//std::cout << "\t\t\t\tBF Method:\n\n";
	//std::cout << "MiddleLength = " << MiddleValueBF << "\tSigma = " << SigmaBF
	//			<<"\tDeltaQ = "<<deltaQBF << std::endl
	//			<< "Confedience interval = {" << intervalBF.first << ", " << intervalBF.second << "}";

	//std::cout << std::endl << std::endl << "\t\t\t\tRS Method:\n\n";
	//std::cout << "MiddleLength = " << MiddleValueRS << "\tSigma = " << SigmaRS
	//	<< "\tDeltaQ = " << deltaQRS << std::endl
	//	<< "Confedience interval = {" << intervalRS.first << ", " << intervalRS.second << "}";

	//std::cout << std::endl << std::endl << "\t\t\t\tG Method:\n\n";
	//std::cout << "MiddleLength = " << MiddleValueG << "\tSigma = " << SigmaG
	//	<< "\tDeltaQ = " << deltaQG << std::endl
	//	<< "Confedience interval = {" << intervalG.first << ", " << intervalG.second << "}";
}


void GamiltonPath::writeToFile(const std::vector<int>& MinLengths, const std::string & name)
{
	std::string title(name);
	std::ofstream statistic;
	statistic.open(title);

	if (!statistic.is_open())
	{
		std::cout << "Error with opening\n";
		exit(1);
	}
	//std::cout << "File open\n";

	for (auto i:MinLengths)
	{
		statistic << i << "\n";
		//std::cout << "\nDone!";
	}
	
	statistic.close();
}

void GamiltonPath::DoMinLenghtVectorBF()
{
	int min = 999; //minimal cost
	int StartVertice = 0;

	do
	{
		if (Vertices[0] != StartVertice)
		{
			std::next_permutation(Vertices.begin(), Vertices.end());
		}

		int sum = 0;
		sum = GetCost(Vertices);
		/*for (size_t i = 0; i < Vertices.size() - 1; i++)
		{
			sum += WeightGraph[Vertices[i]][Vertices[i + 1]];
		}*/

		if (sum<min)
		{
			min = sum;
		}


	} while (std::next_permutation(Vertices.begin(), Vertices.end()));

	MinLengths.push_back(min);
}

void GamiltonPath::DoDeltaQ(float & deltaQ, const float & Sigma, const std::vector<int> & MinLengths)
{
	const float quantil = 1.96; //for 0.95 confidence probability 
	deltaQ = quantil*Sigma / sqrt(MinLengths.size());
}

void GamiltonPath::DoRandVec()
{
	srand(time(NULL));

	std::random_device gen;
	std::mt19937(rd());
	
	int RandSize = 7/*rand() % 4 + 5*/; //number of vertices
	
	//if (debugProperty < RandSize) debugProperty = RandSize; //debug temp

	Vertices.resize(RandSize);     //vector vertices

	for (size_t i = 0; i < Vertices.size(); i++)
	{
		Vertices[i] = i;
	}

	WeightGraph.reserve(RandSize);
	WeightGraph.resize(RandSize);


	//filling graph
	for (size_t i = 0; i < RandSize; i++)
	{
		WeightGraph[i].resize(RandSize);
		for (size_t j = 0; j < RandSize; j++)
		{
			i == j ? WeightGraph[i][j] = 0 : WeightGraph[i][j] = gen() % 16 + 1;
		}
	} 

	//simmetry
	for (size_t i = 0; i < WeightGraph.size(); i++)
	{
		for (size_t j = 0; j < WeightGraph.size(); j++)
		{
			WeightGraph[j][i] = WeightGraph[i][j];
		}
	}

	//general WeightGraph for all
	VerticesSA = VerticesG = VerticesRS = Vertices;
}

void GamiltonPath::DoMinLenghtVectorG()
{
	int sum = 0;

	/*for (size_t i = 0; i < VerticesRS.size() - 1; i++)
	{
		sum += WeightGraph[VerticesG[i]][VerticesG[i + 1]];
	}*/
	sum = GetCost(VerticesG);

	MinLengthsG.push_back(sum);
}

void GamiltonPath::DoMinLenghtVectorRS()
{

	int min = 999; //minimal cost(start)
	int StartVertice = 0;
	int countTimes = VerticesRS.size()*2; //can change
	
	while(countTimes--)
	{
		std::random_shuffle(VerticesRS.begin(), VerticesRS.end());
		if (VerticesRS[0] != StartVertice)
		{
			std::random_shuffle(VerticesRS.begin(), VerticesRS.end());
		}

		int sum = 0;
		for (size_t i = 0; i < VerticesRS.size() - 1; i++)
		{
			sum += WeightGraph[VerticesRS[i]][VerticesRS[i+1]];
		}

		if (sum<min)
		{
			min = sum;	
		}
		

	} 

	MinLengthsRS.push_back(min);
}

void GamiltonPath::DoMinLengthVectorSA()
{
	srand(time(0));
	std::random_device gen;
	std::mt19937(rd());
	
	std::sort(VerticesSA.begin(), VerticesSA.end());

	int costBest = 100;
	int costCurr = 999;
	int costNew = 0;

	float tmpT = temperature;
	float randK = 0;
	float p = 0;
	
	std::vector<int> buffer(2,0);
	
	costCurr = GetCost(VerticesSA);

	while(temperature > 0.01 )
	{
		int randPosFrom = gen() % (VerticesSA.size() - 3) + 1; //get rand position to copy from
		auto FromIt = VerticesSA.begin() + randPosFrom; //make iterator by randPosFrom

		std::copy_n(FromIt, 2, buffer.begin()); //copy a part of vertices in the buffer
		VerticesSA.erase(FromIt, FromIt + 2); //erase elements we coped

		int randPosTo = gen() % (VerticesSA.size() - 1) + 1; //get rand position to insert
		if (abs(randPosFrom - randPosTo) < 2)
			randPosTo += gen() % 1 + 1;

		auto ToIt = VerticesSA.begin() + randPosTo; //make iterator by randPosTo
		auto a = VerticesSA.insert(ToIt, buffer.begin(), buffer.end()); //insert from buffer{} 

		costNew = GetCost(VerticesSA);

		if (costNew - costCurr > 0)
		{
			p = std::exp((costCurr - costNew) / temperature);
			randK = float(rand()%100/100.0);

			if (randK < p)
				costCurr = costNew;
		}
		else if (costNew < costCurr)
		{
			costCurr = costNew;
		}


		if (costCurr < costBest)
			costBest = costCurr;

		temperature *= alpha;
	}
	temperature = tmpT;

	MinLengthsSA.push_back(costBest);
}

int GamiltonPath::GetCost(const std::vector<int>& Vertices)
{
	int sum = 0;

	for (size_t i = 0; i < Vertices.size() - 1; i++)
	{
		sum += WeightGraph[Vertices[i]][Vertices[i + 1]];
	}
	return sum;
}

void GamiltonPath::DoMiddleValue(float & MiddleValue, const std::vector<int> & MinLengths)
{
	MiddleValue = std::accumulate(MinLengths.begin(), MinLengths.end(), 0.0) / MinLengths.size();
}

void GamiltonPath::DoSigma(float & Sigma, const float & MiddleValue, const std::vector<int> & MinLengths)
{
	for (int i : MinLengths)
	{
		Sigma += (i - MiddleValue)*(i - MiddleValue);
	}

	Sigma = sqrt(Sigma / MinLengths.size());
}


void GamiltonPath::DoInterval(std::pair<float,float> & interval, const float & deltaQ, const float & MiddleValue)
{
	interval = { MiddleValue - deltaQ, deltaQ + MiddleValue };
}

