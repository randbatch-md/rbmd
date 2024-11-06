#include <torch/torch.h>
#include <torch/script.h>
#include <algorithm>
#include <iostream>
#include <array>
#include <mace/maceload.h>
#include <torch/jit.h>
#include <chrono>

maceload macetest;

maceload::maceload(){};
void maceload::init(int nodes_, bool v_g, std::string path_, std::string device_name)
{
  torch::jit::getProfilingMode() = false;
  torch::jit::getExecutorMode() = true;
  n_nodes = nodes_;
  std::cout << n_nodes<<" atoms " << std::endl;
  vflag_global = bool(v_g);
  //Kokkos::initialize();
  if (device_name != "cuda" && device_name != "dcu")
  {
    device = c10::Device(torch::kCPU);
  //  using DeviceType = Kokkos::DefaultHostExecutionSpace;
  }
      /*
    if (!torch::cuda::is_available()) {
        std::cout << "CUDA unavailable, setting device type to CPU" << std::endl;
        device = c10::Device(torch::kCPU);
    }
    else {
        device = torch::Device(torch::DeviceType::CUDA);
    }*/
    try {
        model = torch::jit::load(path_, device);
        std::cout << "Loading MACE model from " << path_ << "\" ..."<< std::endl;
    }
    catch (const c10::Error& e) {
        std::cerr << "Error load mace model\n";
    }
    mace_r_max = model.attr("r_max").toTensor().item<double>();
    //double num_inter = model.attr("num_interactions").toTensor().item<double>();
    //std::cout << num_inter <<std::endl;
    auto mace_atom_table = model.attr("atomic_numbers").toTensor();
    n_node_feats = mace_atom_table.numel();
    std::string transsymbol;
    int transnum;
    float transmass;
    for (int a_n = 0; a_n < n_node_feats; ++a_n) {
        transnum = mace_atom_table[a_n].item<int>() -1;
        transmass = mass_periodic_table[transnum];
        mace_feats_mass_table.push_back(transmass);
        transsymbol = periodic_table[transnum];
        mace_feats_table.push_back(transsymbol);
        //std::cout<< transsymbol << std::endl;
        //std::cout << transmass << std::endl;
    }    
    energy = energy.to(device);
    //torch::Tensor positions;
    //torch::Tensor batch;
    //torch::Tensor forces;
    //torch::Tensor node_attrs;
    //torch::Tensor edge_index;
    //torch::Tensor unit_shifts;
    //torch::Tensor shifts;
    //torch::Tensor mask;
}
void maceload::loadatoms(std::vector<std::string> symbollist) {
    /* no change between two timestep */
    batch = torch::zeros({ n_nodes }, torch::dtype(torch::kInt64).device(device));
    forces = torch::empty({ n_nodes, 3 }, torch_float_dtype.device(device));

    ptr[1] = n_nodes;
    ptr = ptr.to(device);
    weight[0] = 1.0;
    weight = weight.to(device);
    node_attrs = torch::zeros({ n_nodes, n_node_feats }, torch_float_dtype.device(device));

    for (int n_a_i = 0; n_a_i < n_nodes; ++n_a_i) {
        //std::string as = symbollist[n_a_i];
        std::string as = symbollist[n_a_i];
        //大小写问题
        try {
            auto iter = std::find(mace_feats_table.begin(), mace_feats_table.end(), as);
            node_attrs[n_a_i][std::distance(mace_feats_table.begin(), iter)] = 1.0;
        }
        catch (const c10::Error& e) {
            std::cerr << "Unsupported element symbol!\n";
        }
    }
    mask = torch::ones({ n_nodes }, torch::dtype(torch::kBool));
   // for (int i = 0; i < n_nodes; ++i) {
   //     mask[i] = true;
   // }
}
void maceload::loadmass(std::vector<float> masslist)
{
  /* no change between two timestep */
  batch = torch::zeros({ n_nodes }, torch::dtype(torch::kInt64).device(device));
  forces = torch::empty({ n_nodes, 3 }, torch_float_dtype.device(device));

  ptr[1] = n_nodes;
  ptr = ptr.to(device);
  weight[0] = 1.0;
  weight = weight.to(device);
  node_attrs = torch::zeros({ n_nodes, n_node_feats }, torch_float_dtype.device(device));
  for (int n_a_i = 0; n_a_i < n_nodes; ++n_a_i)
  {
    float am = masslist[n_a_i];
    float dvalue = 0;
    int dvalueflag = 0;
    for (int n_f_i = 0; n_f_i < n_node_feats; ++n_f_i)
    {
      dvalue = am - mace_feats_mass_table[n_f_i];
      if (dvalue < 0.05 && dvalue > -0.05)
      {
        node_attrs[n_a_i][n_f_i] = 1.0;
        dvalueflag++;
      }
    }
    if (dvalueflag != 1)
    {
      std::cout << "No." << n_a_i
                << " atom have an inaccurate atomic mass.\n Please check or enter atomic symbol."
                << std::endl;
    }
  }
  mask = torch::ones({ n_nodes }, torch::dtype(torch::kBool));
  //for (int i = 0; i < n_nodes; ++i)
  //{
  //  mask[i] = true;
  //}
  //std::cout << "masses" << std::endl;
}
void maceload::loadcell(const Vec3f& box)
{
    cell[0][0] = box[0];
    cell[1][1] = box[1];
    cell[2][2] = box[2];
    /*
    h[0] = xprd;
    h[1] = yprd;
    h[2] = zprd;
    h_inv[0] = 1.0/h[0];
    h_inv[1] = 1.0/h[1];
    h_inv[2] = 1.0/h[2];
    h[3] = yz;
    h[4] = xz;
    h[5] = xy;
    cell[2][1] = h[3];
    cell[2][0] = h[4];
    cell[1][0] = h[5];
    h_inv[3] = -h[3] / (h[1]*h[2]);
    h_inv[4] = (h[3]*h[5] - h[1]*h[4]) / (h[0]*h[1]*h[2]);
    h_inv[5] = -h[5] / (h[0]*h[1]);*/
    cell = cell.to(device);
}
void maceload::loadpositions(const vtkm::cont::ArrayHandle<Vec3f>& position_rbmd)
{
  positions = torch::empty({ n_nodes, 3 }, torch_float_dtype);
    auto read_position = position_rbmd.ReadPortal();
    
    for (size_t i = 0; i < n_nodes; i++)
    {
    positions[i][0] = read_position.Get(i)[0];
    positions[i][1] = read_position.Get(i)[1];
    positions[i][2] = read_position.Get(i)[2];
    }
    //positions = positions.to(device);
    //std::cout << "positions" << std::endl;
}

/*void maceload::loadedges_index(const std::vector<Id> edge0,
                               const std::vector<Id> edge1,
                               const std::vector<Vec3f> _unit_shifts,
                               const std::vector<Vec3f> _shifts)
{
    int64_t index = edge0.size();
    
    // 创建 Kokkos 视图
    auto k_edge_index = Kokkos::View<int64_t**, Kokkos::LayoutRight, DeviceType>("k_edge_index", 2, index);
    auto k_unit_shifts = Kokkos::View<double**, Kokkos::LayoutRight, DeviceType>("k_unit_shifts", index,3);
    auto k_shifts = Kokkos::View<double**, Kokkos::LayoutRight, DeviceType>("k_shifts", index,3);

    // 使用 Kokkos 并行处理
    Kokkos::parallel_for("LoadEdgesIndex", index, KOKKOS_LAMBDA(int i) {
        //k_edge_index(0, i) = edge0[i];
        //k_edge_index(1, i) = edge1[i];
        k_unit_shifts(i, 0) = _unit_shifts[i][0];
        k_unit_shifts(i, 1) = _unit_shifts[i][1];
        k_unit_shifts(i, 2) = _unit_shifts[i][2];
        k_shifts(i, 0) = _shifts[i][0];
        k_shifts(i, 1) = _shifts[i][1];
        k_shifts(i, 2) = _shifts[i][2];
    });

    // 导入 PyTorch 张量
    edge_index = torch::from_blob(
        k_edge_index.data(),
        {2, index},
        torch::TensorOptions().dtype(torch::kInt64).device(device));
    
    unit_shifts = torch::from_blob(
        k_unit_shifts.data(),
        {index, 3},
        torch::TensorOptions().dtype(torch::kF64).device(device));
    
    shifts = torch::from_blob(
        k_shifts.data(),
        {index, 3},
        torch::TensorOptions().dtype(torch::kF64).device(device));
}*/
void maceload::loadedges_index(std::vector<Id> edge0,
                               std::vector<Id> edge1,
                               std::vector<Vec3f> _unit_shifts,
                               std::vector<Vec3f> _shifts)
{
    //edge_index = torch::empty({ 2,n_edges }, torch::dtype(torch::kInt64));
    //unit_shifts = torch::zeros({ n_edges,3 }, torch_float_dtype);
    //shifts = torch::zeros({ n_edges,3 }, torch_float_dtype);
   //float _cut_off = mace_r_max;
    /*
    std::vector<int64_t> edge0;
    std::vector<int64_t> edge1;
    std::vector<double> unit_shifts0;
    std::vector<double> unit_shifts1;
    std::vector<double> unit_shifts2;
    std::vector<double> shifts0;
    std::vector<double> shifts1;
    std::vector<double> shifts2;
    const Real _cut_off_2 = _cut_off * _cut_off;
    for (int ii=0; ii<n_nodes; ++ii) {
        auto p_i = _locator.GetPtsPosition(ii);
        vtkm::Id3 p_i_cell = _locator.InWhichCell(p_i);
        const auto num_cycles = _locator.GetNumCycles();
        for (Id i = -num_cycles; i <= num_cycles; i++)
        {
            for (Id j = -num_cycles; j <= num_cycles; j++)
            {
                for (Id k = -num_cycles; k <= num_cycles; k++)
                {
                    auto neighbor_ijk = p_i_cell + Id3{ i, j, k };
                    auto ijk = _locator.PeriodicIndexOffset(neighbor_ijk);
                    auto coord_offset = _locator.PeriodicCoordOffset(ijk - neighbor_ijk);
                    auto num_pts = _locator.NumberPtsInCell(ijk);

                    for (Id p = 0; p < num_pts; p++)
                    {
                        auto pts_id_j = _locator.PtsInCell(ijk, p);//获取第p个原子id
                        auto p_j = _locator.GetPtsPosition(pts_id_j) - coord_offset;//已经是向量了
                        auto r_ij = p_j - p_i;
                        const Real dis_2 = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
                        if (_cut_off_2 - dis_2 > 0.0001 && dis_2 > 0.0001)
                        {
                            edge0.push_back(ii);
                            edge1.push_back(pts_id_j);                         
                            double shiftxs = std::round(-coord_offset[0] /cell[0][0]);//h012是xyz的长度
                            double shiftys = std::round(-coord_offset[1] / cell[0][0]);
                            double shiftzs = std::round(-coord_offset[2] / cell[0][0]);
                            unit_shifts0.push_back(shiftxs);
                            unit_shifts1.push_back(shiftys);
                            unit_shifts2.push_back(shiftzs);
                            shifts0.push_back(-coord_offset[0]);
                            shifts1.push_back(-coord_offset[1]);
                            shifts2.push_back(-coord_offset[2]);
                        }
                    }
                }
            }
        }
    }*/
    //std::cout << "edges_start" << std::endl;
    //auto start_time1 = std::chrono::high_resolution_clock::now();
    int64_t index = edge0.size();
    //std::cout << index << std::endl;
    if (std::is_same<Id, vtkm::Int64>::value)
    {
      edge_index =
        torch::cat({ torch::from_blob(edge0.data(), { 1, index }, optsint64),
                     torch::from_blob(edge1.data(), { 1, index }, optsint64) },
                   0);
    }
    else if (std::is_same<Id, vtkm::Int32>::value)
    {
      edge_index =
        torch::cat({ torch::from_blob(edge0.data(), { 1, index }, optsint32).contiguous(),
                     torch::from_blob(edge1.data(), { 1, index }, optsint32).contiguous() },
                   0)
          .to(torch::kInt64);
    }
    if (std::is_same<FloatDefault, vtkm::Float64>::value)
    {
      shifts = torch::from_blob(_shifts.data(), { index, 3 }, optsfloat64);
      unit_shifts = torch::from_blob(_unit_shifts.data(), { index, 3 }, optsfloat64);
    }
    else if (std::is_same<FloatDefault, vtkm::Float32>::value)
    {
      shifts = torch::from_blob(_shifts.data(), { index, 3 }, optsfloat32).to(torch::kF64);
      unit_shifts = torch::from_blob(_unit_shifts.data(), { index, 3 }, optsfloat32).to(torch::kF64);
    }
    //edge_index = torch::cat({ torch::from_blob(edge0.data(), {1,index },optsint32).clone(), torch::from_blob(edge1.data(), { 1,index }, optsint32).clone() }, 0).to(torch::kInt64);
    /*edge_index = torch::empty({ 2, index }, torch::dtype(torch::kInt64));
    unit_shifts = torch::zeros({ index, 3 }, torch_float_dtype);
    shifts = torch::zeros({ index, 3 }, torch_float_dtype);

    for (size_t i = 0; i < index; i++)
    {
      edge_index[0][i] =edge0[i] ;
      edge_index[1][i] = edge1[i];
      unit_shifts[i][0] = _unit_shifts[i][0];
      unit_shifts[i][1] = _unit_shifts[i][1];
      unit_shifts[i][2] = _unit_shifts[i][2];
      shifts[i][0] = _shifts[i][0];
      shifts[i][1] = _shifts[i][1];
      shifts[i][2] = _shifts[i][2];
    }*/
    //auto end_time1 = std::chrono::high_resolution_clock::now();
    //auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(end_time1 - start_time1);
    //std::cout << "Time elapsed1: " << duration1.count() << " milliseconds" << std::endl;
    //std::cout << "edges_index" << std::endl;
    /*
    torch::TensorOptions opts = torch::TensorOptions().dtype(torch::kInt64);
    int edges = edge0.size();
    torch::Tensor cat_1 = torch::cat({ torch::from_blob(edge0.data(), {1,edges },opts).clone(), torch::from_blob(edge1.data(), { 1,edges }, opts).clone() }, 0);
    std::cout << cat_1 << std::endl;
    torch::Tensor cat_2 = torch::cat({ torch::from_blob(unit_shifts0.data(), {edges },opts).clone(), torch::from_blob(unit_shifts1.data(), {edges }, opts).clone(), torch::from_blob(unit_shifts2.data(), {edges }, opts).clone() }, 1);
    torch::Tensor cat_3 = torch::cat({ torch::from_blob(shifts0.data(), {edges },opts).clone(), torch::from_blob(shifts1.data(), {edges }, opts).clone(), torch::from_blob(shifts2.data(), {edges }, opts).clone() }, 1);
    /*for (int ii = 0; ii<n_nodes; ++ii) {
    int i = list->ilist[ii];
    double xtmp = atom->x[i][0];
    double ytmp = atom->x[i][1];
    double ztmp = atom->x[i][2];
    int *jlist = list->firstneigh[i];
    int jnum = list->numneigh[i];
    int k = first_edge[ii];
	//jnum是邻居数目，随i每次更新
    for (int jj=0; jj<jnum; ++jj) {
      int j = jlist[jj];
      j &= NEIGHMASK;
      double delx = xtmp - atom->x[j][0];
      double dely = ytmp - atom->x[j][1];
      double delz = ztmp - atom->x[j][2];
      double rsq = delx * delx + dely * dely + delz * delz;
      if (rsq < r_max_squared) {
        edge_index[0][k] = i;
        if (domain_decomposition) {
          edge_index[1][k] = j;
		  //如果有鬼原子，直接链接
		  //否则，按map-tag这个将j_local导出，将差变换后归于unit_shifts,shifts则再经过晶格变换后归入
        } else {
          int j_local = atom->map(atom->tag[j]);
          edge_index[1][k] = j_local;
          double shiftx = atom->x[j][0] - atom->x[j_local][0];
          double shifty = atom->x[j][1] - atom->x[j_local][1];
          double shiftz = atom->x[j][2] - atom->x[j_local][2];
          double shiftxs = std::round(domain->h_inv[0]*shiftx + domain->h_inv[5]*shifty + domain->h_inv[4]*shiftz);
          double shiftys = std::round(domain->h_inv[1]*shifty + domain->h_inv[3]*shiftz);
          double shiftzs = std::round(domain->h_inv[2]*shiftz);
          unit_shifts[k][0] = shiftxs;
          unit_shifts[k][1] = shiftys;
          unit_shifts[k][2] = shiftzs;
          shifts[k][0] = domain->h[0]*shiftxs + domain->h[5]*shiftys + domain->h[4]*shiftzs;
          shifts[k][1] = domain->h[1]*shiftys + domain->h[3]*shiftzs;
          shifts[k][2] = domain->h[2]*shiftzs;
        }
        k++;
      }
    }*/
}


double maceload::energyout(c10::impl::GenericDict ouput)
{
  energy = ouput.at("total_energy_local").toTensor().cpu();
  eng_vdwl = energy.item<double>();
//  std::cout << node_energy_ptr<< std::endl;

  return eng_vdwl;
}


c10::impl::GenericDict maceload::forward() {
  /*if (!torch::cuda::is_available())
   {
        std::cout << "CUDA unavailable, setting device type to CPU" << std::endl;
      device = c10::Device(torch::DeviceType::CPU);
    }
    else {
        device = torch::Device(torch::DeviceType::CUDA);
    }*/
    //std::cout << "forward_start"<< std::endl;
    //auto start_time = std::chrono::high_resolution_clock::now();
    c10::Dict<std::string, torch::Tensor> input;
    try
    {
      //batch = batch.to(device);
      //cell = cell.to(device);
      edge_index = edge_index.to(device);
      //energy = energy.to(device);
      //forces = forces.to(device);
      //node_attrs = node_attrs.to(device);
      positions = positions.to(device);
      //ptr = ptr.to(device);
      shifts = shifts.to(device);
      unit_shifts = unit_shifts.to(device);
      //weight = weight.to(device);
    }
    catch (const c10::Error& e)
    {
      std::cerr << "Error in to device\n";
    }
    /* std::cout << batch.sizes() << std::endl;
    std::cout << cell.sizes() << std::endl;
    std::cout << edge_index.sizes() << std::endl;
    std::cout << energy.sizes() << std::endl;
    std::cout << forces.sizes() << std::endl;
    std::cout << node_attrs.sizes() << std::endl;
    std::cout << positions.sizes() << std::endl;
    std::cout << ptr.sizes() << std::endl;
    std::cout << shifts.sizes() << std::endl;
    std::cout << unit_shifts.sizes() << std::endl;
    std::cout << weight.sizes() << std::endl; 
    std::cout << batch.dtype() << std::endl; 
    std::cout << cell.dtype() << std::endl; 
    std::cout << edge_index.dtype() << std::endl; 
    std::cout << energy.dtype() << std::endl; 
    std::cout << forces.dtype() << std::endl; 
    std::cout << node_attrs.dtype() << std::endl; 
    std::cout << positions.dtype() << std::endl; 
    std::cout << ptr.dtype() << std::endl; 
    std::cout << shifts.dtype() << std::endl; 
    std::cout << unit_shifts.dtype() << std::endl; 
    std::cout << weight.dtype() << std::endl; */
    //std::cout << edge_index.sizes() << std::endl;
    input.insert("batch", batch);
    input.insert("cell", cell);
    input.insert("edge_index", edge_index);
    input.insert("energy", energy);
    input.insert("forces", forces);
    input.insert("node_attrs", node_attrs);
    input.insert("positions", positions);
    input.insert("ptr", ptr);
    input.insert("shifts", shifts);
    input.insert("unit_shifts", unit_shifts);
    input.insert("weight", weight);
    //std::cout << "数据载入完成" << std::endl;

    model.eval();
    c10::impl::GenericDict output = model.forward({ input, mask.to(device), bool(vflag_global) }).toGenericDict();
    //auto end_time = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    //std::cout << "Time elapsed: " << duration.count() << " milliseconds" << std::endl;
    //std::cout << "模型前向完成" << std::endl;
    return output;
}

std::vector<Vec3f> maceload::forcesout(c10::impl::GenericDict ouput)
{
    //auto start_timef = std::chrono::high_resolution_clock::now();
    forces = ouput.at("forces").toTensor().cpu();
    /* if (!forces.is_contiguous())
    {
      std::cout << "s" << std::endl;
      forces = forces.contiguous(); // 确保张量在内存中是连续的
    }*/
    //auto start_timef = std::chrono::high_resolution_clock::now();
    std::vector<Vec3f> _force(n_nodes);
    for (size_t i = 0; i < n_nodes; i++)
    {
      _force[i][0] = forces[i][0].item<float>();
      _force[i][1] = forces[i][1].item<float>();
      _force[i][2] = forces[i][2].item<float>();
    }
    //auto end_timef = std::chrono::high_resolution_clock::now();
    //uto durationf = std::chrono::duration_cast<std::chrono::milliseconds>(end_timef - start_timef);
    //std::cout << "force elapsed: " << durationf.count() << " milliseconds" << std::endl;
  //  std::cout << forces << std::endl;
    return _force;
}
torch::Tensor maceload::virialout(c10::impl::GenericDict ouput) {
    torch::Tensor vir = ouput.at("virials").toTensor().cpu();
    torch::Tensor virial = torch::zeros({6}, torch_float_dtype);
    virial[0] += vir[0][0][0].item<double>();
    virial[1] += vir[0][1][1].item<double>();
    virial[2] += vir[0][2][2].item<double>();
    virial[3] += 0.5 * (vir[0][1][0].item<double>() + vir[0][0][1].item<double>());
    virial[4] += 0.5 * (vir[0][2][0].item<double>() + vir[0][0][2].item<double>());
    virial[5] += 0.5 * (vir[0][2][1].item<double>() + vir[0][1][2].item<double>());
    return virial;
}
torch::Tensor maceload::atomenergy(c10::impl::GenericDict ouput) {
    torch::Tensor node_energy = ouput.at("node_energy").toTensor().cpu();
    torch::Tensor eatom = torch::zeros({ n_nodes }, torch_float_dtype);
    for (int i = 0; i < n_nodes; ++i) {
        eatom[i] = node_energy[i].item<double>();
    }
    return eatom;
}

/*
int main()
{
    maceload r;
    std::string path = "D:/workspace/macemodel/mace_agnesi_small.model-lammps.pt";
    int x = 5;
    int y = 8;
    bool a = false;
    r.init(x, a, path);
    std::cout << r.mace_r_max << std::endl;
    std::vector<std::string> symbollist;
    r.loadatoms(symbollist);

   // r.loadcell(); 

   // r.loadpositions();
    //r.loadedges_index();
    c10::impl::GenericDict x1 = r.forward();
    auto x2 = r.forcesout(x1);
    x2[0][1];

    std::cout << x2<< std::endl;
    std::cout << r.atomenergy(x1) << std::endl;
    std::cout << r.energyout(x1) << std::endl;
}
*/
