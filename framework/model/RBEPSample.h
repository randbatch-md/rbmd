#pragma once
//#define TEST_CODE
// Prepare P samples of \vec{k_l} for RBE force
struct RBEPSAMPLE
{
  Real _alpha;
  Real _Length;
  Vec3f _box;
  Id _P;
  bool _RBE_random;

  Real Compute_S() const
  {
    const Vec3f& factor = Compute_H();
    Real factor_3 = factor[0] * factor[1] * factor[2];
    Real S = factor_3 - 1;
    return S;
  }

  Vec3f Compute_H() const
  {
    Vec3f H = 0.0;
    for (Id i = 0; i < 3; ++i)
    {
      const Real factor = -(_alpha * _box[i] * _box[i]);
      for (int m = -10; m <= 10; m++)
      {
        Real expx = m * m * factor;
        H[i] += vtkm::Exp(expx);
      }
      H[i] *= vtkm::Sqrt(-(factor) / vtkm::Pi());
    }

    return H;
  }

  Real MH_Algorithm(Real m, Real mu, const Vec3f& sigma, Id dimension) const
  {
    Real x_wait = FetchSample_1D(mu, sigma[dimension]);
    Real m_wait = Real(vtkm::Round(x_wait));
    Real Prob = (Distribution_P(m_wait, dimension) / Distribution_P(m, dimension)) *
                (Distribution_q(m, dimension) / Distribution_q(m_wait, dimension));
    Prob = vtkm::Min(Prob, Real(1.0));

    if (_RBE_random)
    {
      Real u = RandomValue<Real>(0.0, 1.0); //random?
      if (u <= Prob)
        m = m_wait;
      return m;
    }
    else
    {
      Real u = 0.5;
      if (u <= Prob)
        m = m_wait;
      return m;
    }
  }

  Real Distribution_P(const Real& x, const Id dimension) const
  {
    Real P_m = vtkm::Exp(-vtkm::Pow(2 * vtkm::Pi() * x / _box[dimension], 2) / (4 * _alpha));
    P_m = P_m / Compute_H()[dimension];
    return P_m;
  }

  Real Distribution_q(const Real& x, const Id dimension) const
  {
    Real q_m;
    if (x == 0)
    {
      q_m = vtkm::ERF((1.0 / 2) /
                  (vtkm::Sqrt(_alpha * vtkm::Pow(_box[dimension], 2) / vtkm::Pow(vtkm::Pi(), 2))));
    }
    else
      q_m = (vtkm::ERF(((1.0 / 2) + vtkm::Abs(x)) /
               (vtkm::Sqrt(_alpha * vtkm::Pow(_box[dimension], 2) / vtkm::Pow(vtkm::Pi(), 2)))) -
             vtkm::ERF((vtkm::Abs(x) - (1.0 / 2)) /
               (vtkm::Sqrt(_alpha * vtkm::Pow(_box[dimension], 2) / vtkm::Pow(vtkm::Pi(), 2))))) /
        2;
    return q_m;
  }

  Real FetchSample_1D(const Real& mu,
                      const Real& sigma) const // Fetch 1D sample from Gaussion contribution
  {
    Real U1, U2, epsilon;
    epsilon = 1e-6;
    if (_RBE_random)
    {
      do
      {
        U1 = RandomValue<Real>(0.0, 1.0);
      } while (U1 < epsilon);
      U2 = RandomValue<Real>(0.0, 1.0);
      Vec2f ChooseSample{ 0.0, 0.0 };
      ChooseSample[0] = sigma * vtkm::Sqrt(-2.0 * vtkm::Log(U1)) * cos(2 * vtkm::Pi() * U2) + mu;
      //ChooseSample[1] = sigma * vtkm::Sqrt(-2.0 * vtkm::Log(U1)) * sin(2 * vtkm::Pi() * U2) + mu;
      return ChooseSample[0];
    }
    else
    {
      U1 = 0.5;
      U2 = 0.5;
      Vec2f ChooseSample{ 0.0, 0.0 };
      ChooseSample[0] = sigma * vtkm::Sqrt(-2.0 * vtkm::Log(U1)) * cos(2 * vtkm::Pi() * U2) + mu;
      //ChooseSample[1] = sigma * vtkm::Sqrt(-2.0 * vtkm::Log(U1)) * sin(2 * vtkm::Pi() * U2) + mu;
      return ChooseSample[0];
    }
  }

  vtkm::cont::ArrayHandle<Vec3f> Fetch_P_Sample(const Real& mu, const Vec3f& sigma) const
  {
    Real epsilonx = 1e-6; // precision
    Vec3f X_0;
    do
    {
      X_0 = { Real(vtkm::Round(FetchSample_1D(mu, sigma[0]))),
              Real(vtkm::Round(FetchSample_1D(mu, sigma[1]))),
              Real(vtkm::Round(FetchSample_1D(mu, sigma[2]))) };
    } while (vtkm::Abs(X_0[0]) < epsilonx && vtkm::Abs(X_0[1]) < epsilonx &&
             vtkm::Abs(X_0[2]) < epsilonx);

    vtkm::cont::ArrayHandle<Vec3f> P_Sample;
    P_Sample.Allocate(_P);
    auto writePortal = P_Sample.WritePortal();
    writePortal.Set(0, X_0);

    for (int i = 1; i < _P; i++)
    {
      Vec3f X_1 = { MH_Algorithm(X_0[0], mu, sigma,0),
                    MH_Algorithm(X_0[1], mu, sigma,1),
                    MH_Algorithm(X_0[2], mu, sigma,2) };
      writePortal.Set(i, X_1);
      X_0 = X_1;
      if (vtkm::Abs(X_1[0]) < epsilonx && vtkm::Abs(X_1[1]) < epsilonx &&
          vtkm::Abs(X_1[2]) < epsilonx)
      {
        i = i - 1;
        continue;
      }
    }
    return P_Sample;
  }
};