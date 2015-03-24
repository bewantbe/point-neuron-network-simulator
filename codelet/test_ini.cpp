// g++ test_ini.cpp -lboost_program_options
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <iterator>
using namespace std;

int main( int argc, char *argv[ ] )
{
  try {
    int opt;
    po::options_description desc("Allowed options");
    desc.add_options()
      ( "help,h", "produce help message" )
      ( "optimization", po::value< int >(&opt)->default_value(10), "optimization level" )
      ( "verbose", po::value< int >()->implicit_value( 1 ), "enable verbosity (optionally specify level)"  )
      ;

    po::variables_map vm;
    po::store( po::command_line_parser( argc, argv ).options( desc ).run(), vm );
    po::notify( vm );


    if ( vm.count( "help" ) )
    {
      cout << "Usage: options_description [options]\n";
      cout << desc;
      return 0;
    }

    if ( vm.count( "verbose" ) )
    {
      cout << "Verbosity enabled.  Level is " << vm[ "verbose" ].as< int >() << "\n";
    }
  }
  catch( std::exception& e )
  {
    cout << e.what() << "\n";
    return 1;
  }    
  return 0;
}
