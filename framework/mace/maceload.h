#pragma once
#include <torch/torch.h>
#include <torch/script.h>
#include <algorithm>
#include <iostream>
#include <array>
#include <vector>
#include <vtkm/Types.h>
#include"Types.h"


//#include "locator/ContPointLocator.h"

class maceload 
{
    public:
        maceload();
        virtual ~maceload()=default;
        //maceload(){};
        void init(int nodes_, bool x_, std::string path_, std::string device_name);
        void loadpositions(const vtkm::cont::ArrayHandle<Vec3f>& position_rbmd); 
        void loadcell(const Vec3f& box);
        //void loadatoms();
        void loadatoms(std::vector<std::string>);
        void loadmass(std::vector<float>);
        void loadedges_index(std::vector<Id> edge0,
                             std::vector<Id> edge1,
                             std::vector<Vec3f> unit_shifts,
                             std::vector<Vec3f> shifts);
        c10::impl::GenericDict forward();
        double energyout(c10::impl::GenericDict ouput);
        std::vector<Vec3f> forcesout(c10::impl::GenericDict ouput);
        torch::Tensor virialout(c10::impl::GenericDict ouput);
        torch::Tensor atomenergy(c10::impl::GenericDict ouput);
        c10::ScalarType torch_float_dtype = torch::kF64;
        double mace_r_max;
        torch::jit::script::Module model;
        int n_nodes, n_edges;
        int vflag_global;
        //ContPointLocator _locator;
    private:

        c10::Device device = torch::Device(torch::DeviceType::CUDA);
        std::vector<std::string> mace_feats_table;
        std::vector<float> mace_feats_mass_table;
        int64_t n_node_feats;
        torch::Tensor cell = torch::zeros({ 3,3 }, torch_float_dtype);
        torch::Tensor weight = torch::empty({ 1 }, torch_float_dtype);
        torch::Tensor energy = torch::empty({ 1 }, torch_float_dtype);
        torch::Tensor positions;


        torch::Tensor batch;
        torch::Tensor forces;
        torch::Tensor ptr = torch::zeros({ 2 }, torch::dtype(torch::kInt64));
        torch::Tensor node_attrs;
        torch::Tensor edge_index;
        torch::Tensor unit_shifts;
        torch::Tensor shifts;
        torch::Tensor mask;
        const std::array<std::string, 118> periodic_table =
        { "H", "He",
         "Li", "Be",                                                              "B",  "C",  "N",  "O",  "F", "Ne",
         "Na", "Mg",                                                             "Al", "Si",  "P",  "S", "Cl", "Ar",
         "K",  "Ca", "Sc", "Ti",  "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
         "Rb", "Sr",  "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te",  "I", "Xe",
         "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
                           "Hf", "Ta",  "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
         "Fr", "Ra", "Ac", "Th", "Pa",  "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
                           "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og" };
        const std::array<float, 118> mass_periodic_table =
        { 1.008, 4.0026,
         6.94  , 9.0122,                                                                                  10.81,  12.011, 14.007,  15.999,  18.998, 20.180,
         22.990, 24.305,                                                                                  26.982, 28.085, 30.974,  32.06,   35.45,  39.948,
         39.098, 40.078, 44.956, 47.867,  50.942, 51.996,  54.938, 55.845, 58.933, 58.693, 63.546, 65.38, 69.723, 72.630, 74.922,  78.971,  79.904, 83.798,
         85.468, 87.62,  88.906, 91.224,  92.906, 95.95,   98,     101.07, 102.91, 106.42, 107.87, 112.41, 114.82, 118.71, 121.76, 127.60,  126.90, 131.29,
         132.91, 137.33, 138.91, 140.12,  140.91, 144.24,  145,    150.36, 151.96, 157.25, 158.93, 162.50, 164.93, 167.26, 168.93, 173.05,  174.97,
                                 178.49, 180.95,  183.84,  186.21, 190.23, 192.22, 195.08, 196.97, 200.59, 204.38, 207.2,  208.98,  209,    210,     222,
         223,    226,    227,    232.04,  231.04, 238.03,  237,    244,    243,    247,    247,    251,    252,    257,    258,     259,    266,
                                 267,    268,     269,     270,    277,    278,    281,    282,    285,    286,    289,    290,     293,    294,     294 };
};


extern maceload macetest;