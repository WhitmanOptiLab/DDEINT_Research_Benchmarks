#ifndef _PLOT_HPP_
#define _PLOT_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

    class Plot
    {
        private:
            void generateGnuplotScript();

            std::string outputfilename;
            std::string title;
            std::string xlabel;
            std::string ylabel;
            std::vector<double> xrange;
            std::vector<double> yrange;
            // (inputfilename, plotType, color, label)
            std::vector<std::tuple<std::string, std::string, std::string, std::string>> plots;

        public:
            Plot(const std::string &outputfilename);
            ~Plot();
            void plot(const std::vector<double> &x, const std::vector<double> &y, const std::string &inputfilename, const std::string &plotType, const std::string &color, const std::string &label);
            void setXrange(const double &min, const double &max);
            void setYrange(const double &min, const double &max);
            void setXlabel(const std::string &xlabel);
            void setYlabel(const std::string &ylabel);
            void setTitle(const std::string &title);
            void show();

            const std::string &getFilename() const;
            const std::string &getTitle() const;
            const std::string &getXlabel() const;
            const std::string &getYlabel() const;
            const std::vector<double> &getXrange() const;
            const std::vector<double> &getYrange() const;
    };


#endif // _PLOT_HPP_