/*
 * overlaps.h
 *
 *  Created on: Jan 24, 2017
 *      Author: isovic
 */

#ifndef SRC_REFACTOR_OVERLAPS_H_
#define SRC_REFACTOR_OVERLAPS_H_

#include <map>
#include <string>
#include <vector>
#include "types.h"

namespace is {

class Overlap {
 public:
  Overlap();

  Overlap(const std::string& line, const OverlapFormat& of, const MapId &q_ids, const MapId &t_ids);

  Overlap(int64_t Aid, int64_t Bid, const std::string &Aname,
          const std::string &Bname, double perc_err, int64_t shared_minmers,
          int64_t Arev, int64_t Astart, int64_t Aend, int64_t Alen,
          int64_t Brev, int64_t Bstart, int64_t Bend, int64_t Blen);
  Overlap(const Overlap& op);
  ~Overlap();
  Overlap& operator=(const Overlap& op);
  void swap(Overlap& op);

  /** Checks the length of the overlap in query and target coordinates.
   * Returns 0 if the ratio of those is within the specified error_rate.
   */
  int CheckConstraints(float error_rate);

  int64_t Aid() const {
    return Aid_;
  }
  int64_t Bid() const {
    return Bid_;
  }
  std::string Aname() const {
    return Aname_;
  }
  std::string Nname() const {
    return Bname_;
  }
  double perc_err() const {
    return perc_err_;
  }
  int64_t shared_minmers() const {
    return shared_minmers_;
  }
  int64_t Arev() const {
    return Arev_;
  }
  int64_t Astart() const {
    return Astart_;
  }
  int64_t Aend() const {
    return Aend_;
  }
  int64_t Alen() const {
    return Alen_;
  }
  int64_t Brev() const {
    return Brev_;
  }
  int64_t Bstart() const {
    return Bstart_;
  }
  int64_t Bend() const {
    return Bend_;
  }
  int64_t Blen() const {
    return Blen_;
  }

 private:
  int ParseMhap_(const std::string &line);
  int ParsePaf_(const std::string &line, const std::map<std::string, int64_t> &q_ids, const std::map<std::string, int64_t> &t_ids);

  int64_t Aid_, Bid_;
  std::string Aname_, Bname_;
  double perc_err_;
  int64_t shared_minmers_;
  int64_t Arev_, Astart_, Aend_, Alen_;  // start is zero-based, end points to a position right after the last inclusive base.
  int64_t Brev_, Bstart_, Bend_, Blen_;
};

class Overlaps {
 public:
  Overlaps(const std::string& path, const OverlapFormat &of, const MapId &q_ids,
           const MapId &t_ids, float error_rate, bool filter_unique);
  ~Overlaps();

  void SortByTargetId();

  const std::vector<Overlap>& overlaps() const {
    return overlaps_;
  }

 private:
  void ParseUnique_(const std::string& path, const OverlapFormat &of,
              const MapId &q_ids, const MapId &t_ids, float error_rate);
  void ParseAll_(const std::string& path, const OverlapFormat &of,
              const MapId &q_ids, const MapId &t_ids, float error_rate);

  std::vector<Overlap> overlaps_;

};

} /* namespace is */

#endif /* SRC_REFACTOR_OVERLAPS_H_ */
