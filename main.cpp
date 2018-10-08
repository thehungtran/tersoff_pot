#include<iostream>
#include<cstdlib>
#include<fstream>
#include<string>
#include<vector>
#include<math.h>
#include<time.h>
#include<sstream>
//#include "pot_en_atom.h"
#include "force_field_tersoff.h"

using namespace std;

int main (int argc, char* argv[])
{
  // prepare system, nearest neighbor

  //ForceFieldAtomistic test1;
  ForceFieldTersoff test1;

  test1.read_potential_file("C.tersoff.params");
  //test1.ini_system("dump_sample.lammpstrj");
  //for (int i = 0; i < 2; i++){
  //ofstream filetype;
  //filetype.open("result.txt");
  //test1.ini_system_bin("dump_sample.lammpstrj112.bin", 0);
  for (int i = 0; i < 5; i++){
  test1.ini_system("dump_sample.lammpstrj", i);}
  //cout << test1.getvalue();
//  test1.compute_force();
//  test1.compute_energy();
  //  filetype.close();
  //}
return 0;
}
