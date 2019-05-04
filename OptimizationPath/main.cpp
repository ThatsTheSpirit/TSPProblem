#include <iostream>
#include "GamiltonPath.h"

using std::cout;
using std::endl;

int main()
{
	GamiltonPath h;
	
	short countTests = 0;
	cout << "Enter the number of tests: ";
	std::cin >> countTests;

	auto begin = std::chrono::steady_clock::now();

	h.SetTemperature(100.0);
	h.SetAlpha(0.99);
	h.DoTest(countTests);
	h.ShowResults();
	
	auto end = std::chrono::steady_clock::now();
	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
	std::cout << "\n\nThe time: " << elapsed_ms.count() << " mS\n";


	cout << endl << endl;
	std::cin.get();
	std::cin.get();


	return 0;
}
