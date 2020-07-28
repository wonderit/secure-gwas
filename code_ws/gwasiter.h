#ifndef __GWASITER_H_
#define __GWASITER_H_

#include "mpc.h"
#include "assert.h"
#include <vector>
#include <NTL/mat_ZZ_p.h>

using namespace std;

class GwasIterator {
public:

//  bool pheno_flag = (mode == GwasIterator::GP_CODE) || (mode == GwasIterator::GMP_CODE);
//  bool missing_flag = (mode == GwasIterator::GM_CODE) || (mode == GwasIterator::GMP_CODE);
  static const int TERM_CODE = 0;
  static const int G_CODE = 1;
//  static const int GP_CODE = 2;
//  static const int GM_CODE = 3;
//  static const int GMP_CODE = 4;

  explicit GwasIterator(MPCEnv& mpc, int pid) {
    this->mpc = &mpc;
    this->pid = pid;
  }

//  void Init(bool pheno_flag, bool missing_flag) {
//    this->pheno_flag = pheno_flag;
//    this->missing_flag = missing_flag;
//
//    index = 0;
//    num_left = Param::NUM_INDS;
//    if (pid == 1) {
//      mpc->SendInt(TransferMode(), 3);
//    }
//
//    cout << "Initialized MIterator" << endl;
//  }

  void Init() {
//    index = 0;
//    num_left = Param::NUM_INDS;
    if (pid == 1) {
      mpc->SendInt(TransferMode(), 3);
    }

    cout << "Initialized Matrix Multiply Iterator" << endl;
  }

  void Terminate() {
    if (pid == 1) {
      mpc->SendInt(TERM_CODE, 3);
    }
  }

  int TransferMode() {
    return G_CODE;
  }

  bool NotDone() {
    return num_left > 0;
  }

  /* Return genotypes in dosage format (0, 1, or 2), assume no missing */
//  void GetNextGP(Vec<ZZ_p>& g, Vec<ZZ_p>& p) {
//    assert(TransferMode() == GP_CODE);
//    Vec<ZZ_p> m;
//    Mat<ZZ_p> gmat;
//    GetNextAux(gmat, m, p);
//    g = gmat[0];
//  }
  void GetNextG(Vec<ZZ_p>& g) {
    assert(TransferMode() == G_CODE);
    Vec<ZZ_p> p;
    Vec<ZZ_p> m;
    Mat<ZZ_p> gmat;
    GetNextAux(gmat, m, p);
    g = gmat[0];
  }

  void GetM1M2(Mat<ZZ_p>& m1, Mat<ZZ_p>& m2) {
    assert(TransferMode() == G_CODE);
    if (pid == 2) {
      mpc->ReceiveMat(m1, 3, Param::M1_NUM_ROW, Param::M1_NUM_COL);
      mpc->ReceiveMat(m2, 3, Param::M2_NUM_ROW, Param::M2_NUM_COL);
    } else if (pid == 1) {
      mpc->SwitchSeed(3);
      MPCEnv::RandMat(m1, Param::M1_NUM_ROW, Param::M1_NUM_COL);
      MPCEnv::RandMat(m2, Param::M2_NUM_ROW, Param::M2_NUM_COL);
//
//      if (pheno_flag) {
//        MPCEnv::RandVec(p, 1 + Param::NUM_COVS);
//      }
//      if (missing_flag) {
//        MPCEnv::RandMat(g, 3, Param::NUM_SNPS);
//        MPCEnv::RandVec(m, Param::NUM_SNPS);
//      } else {
//        MPCEnv::RandMat(g, 1, Param::NUM_SNPS);
//      }
      mpc->RestoreSeed();
    }
  }

//  /* Return genotype probabilities with missingness information */
//  void GetNextGM(Mat<ZZ_p>& g, Vec<ZZ_p>& m) {
//    assert(TransferMode() == GM_CODE);
//    Vec<ZZ_p> p;
//    GetNextAux(g, m, p);
//  }
//  void GetNextGMP(Mat<ZZ_p>& g, Vec<ZZ_p>& m, Vec<ZZ_p>& p) {
//    assert(TransferMode() == GMP_CODE);
//    GetNextAux(g, m, p);
//  }

private:
  MPCEnv *mpc;
  int pid;
  int num_left;
  bool pheno_flag;
  bool missing_flag;

  void GetNext(Mat<ZZ_p>& m1, Mat<ZZ_p>& m2) {
    assert(NotDone());
    m1.SetDims(Param::M1_NUM_ROW, Param::M1_NUM_COL);
    m2.SetDims(Param::M2_NUM_ROW, Param::M2_NUM_COL);

    if (pid == 2) {
      mpc->ReceiveMat(m1, 3, Param::M1_NUM_ROW, Param::M1_NUM_COL);
      mpc->ReceiveMat(m2, 3, Param::M2_NUM_ROW, Param::M2_NUM_COL);
    } else if (pid == 1) {
      mpc->SwitchSeed(3);
//      MPCEnv::RandMat(g, 3, Param::NUM_SNPS);
      MPCEnv::RandMat(m1, Param::M1_NUM_ROW, Param::M1_NUM_COL);
      MPCEnv::RandMat(m2, Param::M2_NUM_ROW, Param::M2_NUM_COL);
      mpc->RestoreSeed();
    }

//    index++;
//    num_left--;
  }

  void GetNextAux(Mat<ZZ_p>& g, Vec<ZZ_p>& m, Vec<ZZ_p>& p) {
    assert(NotDone());

    if (missing_flag) {
      g.SetDims(3, Param::NUM_SNPS);
      m.SetLength(Param::NUM_SNPS);
    } else {
      g.SetDims(1, Param::NUM_SNPS);
    }

    if (pid == 2) {
      if (pheno_flag) {
        mpc->ReceiveVec(p, 3, 1 + Param::NUM_COVS);
      }
      if (missing_flag) {
        mpc->ReceiveMat(g, 3, 3, Param::NUM_SNPS);
        mpc->ReceiveVec(m, 3, Param::NUM_SNPS);
      } else {
        mpc->ReceiveMat(g, 3, 1, Param::NUM_SNPS);
      }
    } else if (pid == 1) {
      mpc->SwitchSeed(3);
      if (pheno_flag) {
        MPCEnv::RandVec(p, 1 + Param::NUM_COVS);
      }
      if (missing_flag) {
        MPCEnv::RandMat(g, 3, Param::NUM_SNPS);
        MPCEnv::RandVec(m, Param::NUM_SNPS);
      } else {
        MPCEnv::RandMat(g, 1, Param::NUM_SNPS);
      }
      mpc->RestoreSeed();
    }

//    index++;
//    num_left--;
  }
};

#endif
