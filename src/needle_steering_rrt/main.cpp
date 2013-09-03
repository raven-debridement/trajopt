#include "interface.h"

// Main driver routine for the wrapper program
int main(int argc, char* argv[])
{
	// Initialize planner
	Interface iface(PLANNER3D);
	
	iface.runNeedleSteering();

	// Exit
	int num;
	std::cin >> num;

	return 0;
}
