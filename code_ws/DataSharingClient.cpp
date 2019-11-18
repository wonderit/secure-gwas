
#include "connect.h"
#include "mpc.h"
#include "protocol.h"
#include "util.h"
#include "gwasiter.h"
#include "NTL/ZZ_p.h"

#include <cstdlib>
#include <fstream>
#include <map>
#include <iostream>
#include <sstream>
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"

using namespace NTL;
using namespace std;

bool send_stream(string data_dir, MPCEnv &mpc, int mode) {
  string m1_file = data_dir + "m1.txt";
  string m2_file = data_dir + "m2.txt";

  ifstream fin_m1(m1_file.c_str());
  ifstream fin_m2(m2_file.c_str());

  if (!fin_m1.is_open()) {
    spdlog::error("Error: could not open m1 {}", m1_file);
    return false;
  }

  if (!fin_m2.is_open()) {
    spdlog::error("Error: could not open m2 {}", m2_file);
    return false;
  }


//  string geno_file = data_dir + "geno.txt";
//  string pheno_file = data_dir + "pheno.txt";
//  string cov_file = data_dir + "cov.txt";

//  bool pheno_flag = true;
//  bool missing_flag = false;
//  bool pheno_flag = (mode == GwasIterator::GP_CODE) || (mode == GwasIterator::GMP_CODE);
//  bool missing_flag = (mode == GwasIterator::GM_CODE) || (mode == GwasIterator::GMP_CODE);

//  ifstream fin_geno(geno_file.c_str());
//  ifstream fin_pheno(pheno_file.c_str());
//  ifstream fin_cov(cov_file.c_str());
//
//  if (!fin_geno.is_open()) {
//    cout << "Error: could not open " << geno_file << endl;
//    return false;
//  }
//
//  if (pheno_flag && !fin_pheno.is_open()) {
//    cout << "Error: could not open " << pheno_file << endl;
//    return false;
//  }
//
//  if (!fin_cov.is_open()) {
//    cout << "Error: could not open " << cov_file << endl;
//    return false;
//  }

  long val ,val2;
  string line1;
  string line2;

  uint32_t m1_lineno = 0;
  uint32_t m2_lineno = 0;

  Mat<ZZ_p> m1;
  Mat<ZZ_p> m2;
  Mat<long> lm1;
  Mat<long> lm2;
  lm1.SetDims(Param::M1_NUM_ROW, Param::M1_NUM_COL);
  lm2.SetDims(Param::M2_NUM_ROW, Param::M2_NUM_COL);

//  m1
  while (getline(fin_m1, line1)) {
    istringstream iss_m1(line1);
    for (int j = 0; j < Param::M1_NUM_COL; j++) {
      iss_m1 >> val;
      lm1[m1_lineno][j] = val;
      spdlog::info("M1 ({},{}): {}", m1_lineno, j, val );
    }
    m1_lineno++;
  }

//  m2
  while (getline(fin_m2, line2)) {
    istringstream iss_m2(line2);
    for (int k = 0; k < Param::M2_NUM_COL; k++) {
      iss_m2 >> val2;
      lm2[m2_lineno][k] = val2;
      spdlog::info("M2 ({},{}): {}", m2_lineno, k, val2 );
    }
    m2_lineno++;
  }
  Mat<ZZ_p> rm1;
  Mat<ZZ_p> rm2;
  mpc.SwitchSeed(1);
  MPCEnv::RandMat(rm1, Param::M1_NUM_ROW, Param::M1_NUM_COL);
  MPCEnv::RandMat(rm2, Param::M2_NUM_ROW, Param::M2_NUM_COL);
  mpc.RestoreSeed();

  spdlog::info("print m1 & m2");
  PrintMat(lm1);
  PrintMat(lm2);
  ZZ_p fp_one;
  IntToFP(m1, lm1, Param::NBIT_K, Param::NBIT_F);
  IntToFP(m2, lm2, Param::NBIT_K, Param::NBIT_F);
//  m1 *= fp_one;
//  m2 *= fp_one;

  spdlog::info("print m1 & m2 after fp_one");
  PrintMat(m1);
  PrintMat(m2);

//  start
  m1 -= rm1;
  m2 -= rm2;
  mpc.SendMat(m1, 2);
  mpc.SendMat(m2, 2);

  spdlog::info("print m1 - rm1");
  PrintMat(m1);
  spdlog::info("print m2 - rm2");
  PrintMat(m2);


  spdlog::info("print rm1");
  PrintMat(rm1);
  spdlog::info("print rm2");
  PrintMat(rm2);
//
// end
//  spdlog::info("print m1 + rm1");
//  m1 += rm1;
//  PrintMat(m1);
//  spdlog::info("print m2 + rm2");
//  m2 += rm2;
//  PrintMat(m2);




  // Generate masks

//    Vec<ZZ_p> rm;
//    if (missing_flag) {
//      MPCEnv::RandMat(rg, 3, Param::NUM_SNPS);
//      MPCEnv::RandVec(rm, Param::NUM_SNPS);
//    } else {
//      MPCEnv::RandMat(rg, 1, Param::NUM_SNPS);
//    }
//    mpc.RestoreSeed();
//
//    // Send masked data
//    g -= rg;
//    mpc.SendMat(g, 2);

//  while (getline(fin_m1, line)) {
//    istringstream iss_m1(line);
//
////    if (pheno_flag) {
//    if (true) {
//      Vec<ZZ_p> p;
//      p.SetLength(1 + Param::NUM_COVS);
//
//      fin_pheno >> p[0];
//
////      if (!getline(fin_m2, line)) {
////        cout << cov_file << " has fewer lines than expected" << endl;
////        return false;
////      }
//
//      istringstream iss_m2(line);
//      for (int j = 0; j < Param::NUM_COVS; j++) {
//        iss_cov >> val;
//        p[j + 1] = ZZ_p(val);
//      }
//
//      Vec<ZZ_p> rp;
//      mpc.SwitchSeed(1);
//      MPCEnv::RandVec(rp, 1 + Param::NUM_COVS);
//      mpc.RestoreSeed();
//
//      p -= rp;
//      mpc.SendVec(p, 2);
//    } else {
//      iss_geno >> val;
//    }
//
//    Mat<ZZ_p> g;
//    Vec<ZZ_p> m;
//    if (true) {
////    if (missing_flag) {
//      Init(g, 3, Param::NUM_SNPS);
//      Init(m, Param::NUM_SNPS);
//    } else {
//      Init(g, 1, Param::NUM_SNPS);
//    }
//
//    // Read from file
//    for (int j = 0; j < Param::NUM_SNPS; j++) {
//      string str;
//      iss_geno >> str;
//      if (str == "NA" || str == "-1") {
//        val = -1;
//      } else if (str == "0") {
//        val = 0;
//      } else if (str == "1") {
//        val = 1;
//      } else if (str == "2") {
//        val = 2;
//      } else {
//        cout << "Error: unknown value in dosage matrix (" << str << ")" << endl;
//        return false;
//      }
//      if (true) {
////      if (missing_flag) {
//        if (val == 0) {
//          g[0][j] = 1;
//        } else if (val == 1) {
//          g[1][j] = 1;
//        } else if (val == 2) {
//          g[2][j] = 1;
//        } else {
//          m[j] = 1;
//        }
//      } else {
//        g[0][j] = ZZ_p(val);
//      }
//    }
//
//    // Generate masks
//    Mat<ZZ_p> rg;
//    Vec<ZZ_p> rm;
//    mpc.SwitchSeed(1);
//    if (missing_flag) {
//      MPCEnv::RandMat(rg, 3, Param::NUM_SNPS);
//      MPCEnv::RandVec(rm, Param::NUM_SNPS);
//    } else {
//      MPCEnv::RandMat(rg, 1, Param::NUM_SNPS);
//    }
//    mpc.RestoreSeed();
//
//    // Send masked data
//    g -= rg;
//    mpc.SendMat(g, 2);
//
//    if (true) {
////    if (missing_flag) {
//      m -= rm;
//      mpc.SendVec(m, 2);
//    }
//
//    lineno++;
//  }
//
//  if (lineno != Param::M1_NUM_ROW) {
//    cout << "Error: data matrix does not have M1_NUM_ROW rows" << endl;
//    return false;
//  }

//  fin_geno.close();
//  fin_pheno.close();
//  fin_cov.close();
  fin_m1.close();
  fin_m2.close();

  return true;
}

int main(int argc, char **argv) {
  if (argc < 3) {
    cout << "Usage: DataSharingClient party_id param_file [data_dir (for P3/SP)]" << endl;
    return 1;
  }

  string pid_str(argv[1]);

  spdlog::info(pid_str);

  int pid;
  if (!Param::Convert(pid_str, pid, "party_id") || pid < 0 || pid > 3) {
    cout << "Error: party_id should be 0, 1, 2, or 3" << endl;
    return 1;
  }

  if (!Param::ParseFile(argv[2])) {
    cout << "Could not finish parsing parameter file" << endl;
    return 1;
  }

  string data_dir;
  if (pid == 3) {
    if (argc < 4) {
      cout << "Error: for P3/SP, data directory should be provided as the last argument" << endl;
      return 1;
    }

    data_dir = argv[3];
    if (data_dir[data_dir.size() - 1] != '/') {
      data_dir += "/";
    }

    cout << "Data directory: " << data_dir << endl;
  }

  vector<pair<int, int> > pairs;
  pairs.push_back(make_pair(0, 1));
  pairs.push_back(make_pair(0, 2));
  pairs.push_back(make_pair(1, 2));
  pairs.push_back(make_pair(1, 3));
  pairs.push_back(make_pair(2, 3));

//  spdlog::info("Pairs : {}", pairs.data()->first);

  /* Initialize MPC environment */
  MPCEnv mpc;
  if (!mpc.Initialize(pid, pairs)) {
    cout << "MPC environment initialization failed" << endl;
    return 1;
  }

  bool success;
  if (pid < 3) {
    success = data_sharing_protocol(mpc, pid);
  } else {
    /* Stream data upon request */
    int signal = mpc.ReceiveInt(1);
    spdlog::info("Signal : {}", signal);

    while (signal != GwasIterator::TERM_CODE) {
      success = send_stream(data_dir, mpc, signal);
      if (!success) break;

      signal = mpc.ReceiveInt(1);
      spdlog::info("Signal in While Loop : {}", signal);
    }

    cout << "Done with streaming data" << endl;
    success = true;
  }

  // This is here just to keep P0 online until the end for data transfer
  // In practice, P0 would send data in advance before each phase and go offline
  if (pid == 0) {
    mpc.ReceiveBool(2);
  } else if (pid == 2) {
    mpc.SendBool(true, 0);
  }

  mpc.CleanUp();

  if (success) {
    cout << "Protocol successfully completed" << endl;
    return 0;
  } else {
    cout << "Protocol abnormally terminated" << endl;
    return 1;
  }
}
