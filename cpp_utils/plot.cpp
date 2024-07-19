#include "plot.hpp"

    Plot::Plot(const std::string &outputfilename)
    {
        // system("mkdir -p gnustuff");
        system("mkdir -p plots");
        this->outputfilename = outputfilename;
    }

    Plot::~Plot()
    {
        
    }

    void Plot::plot(const std::vector<double> &x, const std::vector<double> &y, const std::string &inputfilename, const std::string &plotType, const std::string &color, const std::string &label)
    {
        // add to plots
        plots.push_back(std::make_tuple(inputfilename, plotType, color, label));
        // create gnustuff directory if it does not exist
        std::string ouput_dest = inputfilename;
        // write to file
        std::ofstream file(ouput_dest);
        for (size_t i = 0; i < x.size(); ++i)
        {
            file << x[i] << " " << y[i] << "\n";
        }
        file.close();
    }

    void Plot::setXrange(const double &min, const double &max)
    {
        this->xrange = {min, max};
    }
    
    void Plot::setYrange(const double &min, const double &max)
    {
        this->yrange = {min, max};
    }
    
    void Plot::setXlabel(const std::string &xlabel)
    {
        this->xlabel = xlabel;
    }
    
    void Plot::setYlabel(const std::string &ylabel)
    {
        this->ylabel = ylabel;
    }
    
    void Plot::setTitle(const std::string &title)
    {   
        this->title = title;
    }
    
    void Plot::show()
    {
        std::string scriptFilename = outputfilename + ".gp";
        std::ofstream script(scriptFilename);
        // check if file is open
        if (!script.is_open())
        {
            std::cerr << "Error: could not open file " << scriptFilename << std::endl;
            return;
        }
        script << "set terminal pngcairo size 1050,750 enhanced font 'serif,10'\n";
        script << "set output '" << "plots/" << outputfilename << ".png'\n";
        if (!xlabel.empty())
            script << "set xlabel '" << xlabel << "'\n";
        if (!ylabel.empty())
            script << "set ylabel '" << ylabel << "'\n";
        if (!xrange.empty())
            script << "set xrange [" << xrange[0] << ":" << xrange[1] << "]\n";
        if (!yrange.empty())
            script << "set yrange [" << yrange[0] << ":" << yrange[1] << "]\n";
        if (!title.empty())
            script << "set title '" << title << "'\n";
        script << "set grid\n";
        script << "plot ";
        for (size_t i = 0; i < plots.size(); ++i)
        {
            if (std::get<1>(plots[i]) == "lines")
                script << "'" << std::get<0>(plots[i]) << "' with " << std::get<1>(plots[i]) << " linecolor rgb '" << std::get<2>(plots[i]) << "' title '" << std::get<3>(plots[i]) << "'";
            else if (std::get<1>(plots[i]) == "points")
                script << "'" << std::get<0>(plots[i]) << "' with " << std::get<1>(plots[i]) << " pointtype 7 linecolor rgb '" << std::get<2>(plots[i]) << "' title '" << std::get<3>(plots[i]) << "'";
            if (i < plots.size() - 1)
                script << ",\\\n";
        }
        script.close();


        // run gnuplot
        std::string command = "gnuplot " + scriptFilename;
        system(command.c_str());

        // remove the script file
        std::string remove_command = "rm " + scriptFilename;
        system(remove_command.c_str());

        // remove the data files
        for (size_t i = 0; i < plots.size(); ++i)
        {
            std::string remove_command = "rm " + std::get<0>(plots[i]);
            system(remove_command.c_str());
        }
    }
    
    const std::string &Plot::getFilename() const
    {
        return outputfilename;
    }
    
    const std::string &Plot::getTitle() const
    {
        return title;
    }
    
    const std::string &Plot::getXlabel() const
    {
        return xlabel;
    }
    
    const std::string &Plot::getYlabel() const
    {
        return ylabel;
    }
    
    const std::vector<double> &Plot::getXrange() const
    {
        return xrange;
    }
    
    const std::vector<double> &Plot::getYrange() const
    {
        return yrange;
    }
    