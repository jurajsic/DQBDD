#include <easylogging++.hpp>
#include <formula.hpp>

INITIALIZE_EASYLOGGINGPP

int main(int, char **argv)
{
    hqspre::Formula formula;

    // Configure logging
    el::Configurations defaultConf;
    defaultConf.setToDefault();
    defaultConf.setGlobally(el::ConfigurationType::Enabled, "true");
    defaultConf.setGlobally(el::ConfigurationType::Format, "[HQSpre] %msg");
    defaultConf.set(el::Level::Verbose, el::ConfigurationType::Format, "[HQSpre] %msg");
    defaultConf.setGlobally(el::ConfigurationType::ToFile, "false");
    defaultConf.setGlobally(el::ConfigurationType::ToStandardOutput, "true");
    el::Loggers::reconfigureAllLoggers(defaultConf);
    el::Loggers::setVerboseLevel(0);

    formula.settings().consistency_check = false;

    // Parse the file
    std::string in_name(argv[1]);
    std::ifstream in(in_name);
    in >> formula;
    in.close();

    // do the preprocessing magic
    try {
        formula.settings().bla              = false;
        formula.settings().ble              = false;
        formula.settings().pure_sat_timeout = 1000;
        formula.preprocess();
        //formula.printStatistics();
    } catch (hqspre::SATException&) {
        std::cout << "p cnf 0 0\n" << std::endl;
        return 10;
    } catch (hqspre::UNSATException&) {
        std::cout << "p cnf 0 1\n0\n" << std::endl;
        return 20;
    }

    std::cout << formula;
    return 30;
}