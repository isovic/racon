/*
 * overlaps.cc
 *
 *  Created on: Jan 24, 2017
 *      Author: isovic
 */

#include "overlaps.h"
#include <assert.h>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include "types.h"
#include "log_system/log_system.h"

namespace is {

Overlap::Overlap() : Aid_(0), Bid_(0), Aname_(""), Bname_(""), perc_err_(0.0), shared_minmers_(0), Arev_(0), Astart_(0), Aend_(0), Alen_(0), Brev_(0), Bstart_(0), Bend_(0), Blen_(0) {

}

Overlap::Overlap(int64_t Aid, int64_t Bid, const std::string &Aname,
                 const std::string &Bname, double perc_err,
                 int64_t shared_minmers, int64_t Arev, int64_t Astart,
                 int64_t Aend, int64_t Alen, int64_t Brev, int64_t Bstart,
                 int64_t Bend, int64_t Blen) {
  Aid_ = Aid;
  Bid_ = Bid;
  Aname_ = Aname;
  Bname_ = Bname;
  perc_err_ = perc_err;
  shared_minmers_ = shared_minmers;
  Arev_ = Arev;
  Astart_ = Astart;
  Aend_ = Aend;
  Alen_ = Alen;
  Brev_ = Brev;
  Bstart_ = Bstart;
  Bend_ = Bend;
  Blen_ = Blen;
}

Overlap::~Overlap() {
}

Overlap::Overlap(const Overlap& op) {
  Aid_ = op.Aid_;
  Bid_ = op.Bid_;
  Aname_ = op.Aname_;
  Bname_ = op.Bname_;
  perc_err_ = op.perc_err_;
  shared_minmers_ = op.shared_minmers_;
  Arev_ = op.Arev_;
  Astart_ = op.Astart_;
  Aend_ = op.Aend_;
  Alen_ = op.Alen_;
  Brev_ = op.Brev_;
  Bstart_ = op.Bstart_;
  Bend_ = op.Bend_;
  Blen_ = op.Blen_;
}

Overlap& Overlap::operator =(const Overlap& op) {
  if (&op == this) {
    return *this;
  }  // Check if same.
  Overlap temp(op);                   // Thread safe.
  swap(temp);                         // Swaaaaaap!
  return *this;                       // Here you go!
}

void Overlap::swap(Overlap& op) {
  std::swap(this->Aid_, op.Aid_);
  std::swap(this->Bid_, op.Bid_);
  std::swap(this->Aname_, op.Aname_);
  std::swap(this->Bname_, op.Bname_);
  std::swap(this->perc_err_, op.perc_err_);
  std::swap(this->shared_minmers_, op.shared_minmers_);
  std::swap(this->Arev_, op.Arev_);
  std::swap(this->Astart_, op.Astart_);
  std::swap(this->Aend_, op.Aend_);
  std::swap(this->Alen_, op.Alen_);
  std::swap(this->Brev_, op.Brev_);
  std::swap(this->Bstart_, op.Bstart_);
  std::swap(this->Bend_, op.Bend_);
  std::swap(this->Blen_, op.Blen_);
}

Overlap::Overlap(const std::string& line, const OverlapFormat& of,
                 const MapId& q_ids, const MapId& t_ids) {
  assert(of.isPaf() || of.isMhap());
  if (of.isPaf()) {
    ParsePaf_(line, q_ids, t_ids);
  } else if (of.isMhap()) {
    ParseMhap_(line);
  }
}

int Overlap::CheckConstraints(float error_rate) {
  double Adist = Aend_ - Astart_;
  double Bdist = Bend_ - Bstart_;
  double ratio = (Adist > Bdist) ? (1.0f - Bdist / Adist) : (1.0f - Adist/Bdist);
  if (ratio > error_rate) { return 1; }
  return 0;
}

int Overlap::ParseMhap_(const std::string &line) {
  std::istringstream iss(line);
  if (!(iss >> Aid_ >> Bid_ >> perc_err_ >> shared_minmers_ >> Arev_ >> Astart_ >> Aend_ >> Alen_ >> Brev_ >> Bstart_ >> Bend_ >> Blen_)) {
    ERROR_REPORT(ERR_UNEXPECTED_VALUE, "Overlaps are not formatted in the MHAP format. Exiting.");
  }
  Aname_ = std::to_string(Aid_);
  Bname_ = std::to_string(Bid_);
  return 0;
}

int Overlap::ParsePaf_(const std::string &line, const MapId &q_ids, const MapId &t_ids) {
  std::istringstream iss(line);
  std::string tempAstrand, cm;
  int64_t num_residue_matches=0, aln_block_len=0, mapq=0;
  if (!(iss >> Aname_ >> Alen_ >> Astart_ >> Aend_ >> tempAstrand >> Bname_ >> Blen_ >> Bstart_ >> Bend_ >> num_residue_matches >> aln_block_len >> mapq >> cm)) {
    ERROR_REPORT(ERR_UNEXPECTED_VALUE, "Overlaps are not formatted in the PAF format. Exiting.");
  }

  if (tempAstrand == "+") {
    Arev_ = 0; Brev_ = 0;
  } else {
    Arev_ = 0; Brev_ = 1;
  }

  // Convert the scores.
  perc_err_ = ((double) num_residue_matches) / ((double) aln_block_len);
  std::size_t found = cm.find_last_of(":");
  sscanf (cm.substr(found+1).c_str(), "%ld", &shared_minmers_);

  // Find the ID of A.
  auto it_a = q_ids.find(Aname_);
  if (it_a == q_ids.end()) {
    ERROR_REPORT(ERR_UNEXPECTED_VALUE, "Could not find qname '%s' in the input reads file! Exiting.", Aname_.c_str());
  }
  Aid_ = it_a->second + 1;

  // Find the ID of B.
  auto it_b = t_ids.find(Bname_);
  if (it_b == t_ids.end()) {
    ERROR_REPORT(ERR_UNEXPECTED_VALUE, "Could not find qname '%s' in the input contigs file! Exiting.", Bname_.c_str());
  }
  Bid_ = it_b->second + 1;

  return 0;
}



Overlaps::Overlaps(const std::string& path, const OverlapFormat &of,
                   const MapId &q_ids, const MapId &t_ids, float error_rate, bool filter_unique) {
  Parse_(path, of, q_ids, t_ids, error_rate, filter_unique);
}

Overlaps::~Overlaps() {

}

void Overlaps::SortByTargetId() {
  std::sort(overlaps_.begin(), overlaps_.end(),
            [](const Overlap &a, const Overlap &b){ return (a.Bid() < b.Bid()); } );
}

void Overlaps::Parse_(const std::string& path, const OverlapFormat& of,
                      const MapId& q_ids, const MapId& t_ids,
                      float error_rate, bool filter_unique) {
  assert(of.isPaf() || of.isMhap());

  std::unordered_map<int64_t, Overlap> fmap;     // Filtering map.

  std::ifstream infile;
  if (path != "-") {                    // Read from disk.
    infile.open(path.c_str());
    if (!infile) {
      std::cerr << "ERROR: Could not open file '" << path << "' for reading! Exiting." << std::endl;
      exit(1);
    }
  }
  std::istream& input =(path != "-") ? infile : std::cin;       // Choose the input stream.

  std::string line;
  while (std::getline(input, line))
  {
    if (line.size() == 0) { continue; }
    std::istringstream iss(line);

    Overlap overlap(line, of, q_ids, t_ids);      // Parse the line and check if everything went fine.

    if (!overlap.CheckConstraints(error_rate)) {  // Do simple filtering of erroneous overlaps.
      continue;
    }

    // Look-up the A sequence in the map. Update if it's better, or add if
    // nothing yet exists for this A.
    auto it = fmap.find(overlap.Aid());
    if (it == fmap.end() || overlap.shared_minmers() > it->second.shared_minmers()) {
      fmap[overlap.Aid()] = overlap;
    }
  }

  if (infile.is_open()) {
    infile.close();
  }

  overlaps_.clear();
  overlaps_.reserve(fmap.size());
  for (auto it = fmap.begin(); it != fmap.end(); it++) {
    overlaps_.push_back(it->second);
  }
}

} /* namespace is */
