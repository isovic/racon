/*
 * window.h
 *
 *  Created on: Jan 29, 2017
 *      Author: isovic
 */

#ifndef SRC_REFACTOR_WINDOW_H_
#define SRC_REFACTOR_WINDOW_H_

#include <stdint.h>
#include <vector>

#include "types.h"

namespace is {

/** A container object holding info about a single query overlapping
 *  a particular window.
*/
class WindowEntry {
 public:
  WindowEntry() : overlap_id_(-1) { }
  WindowEntry(int64_t overlap_id, int64_t qstart, int64_t qend, int64_t tstart, int64_t tend) :
   			overlap_id_(overlap_id) {
    query_.start = qstart;  query_.end = qend;
    target_.start = tstart; target_.end = tend;
  }
  ~WindowEntry() { }

  int64_t overlap_id() const {
  	return overlap_id_;
  }
  Range query() const {
  	return query_;
  }
  Range target() const {
  	return target_;
  }

 private:
  int64_t overlap_id_ = - 1;
  Range query_, target_;					// Query and target start and end locations.
};

class Window {
 public:
  Window(int64_t target_id);
  ~Window();

  /** Adds info for a part of a sequence overlapping this particular window.
   * @overlap_id ID of the overlap used to generate this info.
   * @qstart Query
   */
  void add(int64_t overlap_id, int64_t qstart, int64_t qend, int64_t tstart, int64_t tend) {
  	entries_.emplace_back(WindowEntry(overlap_id, qstart, qend, tstart, tend));
  }

  const std::vector<WindowEntry>& entries() const {
  	return entries_;
  }

  int64_t target_id() const {
    return target_id_;
  }

  int64_t start() const {
  	return start_;
  }

  void start(int64_t new_start) {
  	start_ = new_start;
  }

  int64_t end() const {
  	return end_;
  }

  void end(int64_t new_end) {
  	end_ = new_end;
  }

 private:
//  std::vector<
  std::vector<WindowEntry> entries_;
  int64_t target_id_;
  int64_t start_;
  int64_t end_;
};

} /* namespace is */

#endif /* SRC_REFACTOR_WINDOW_H_ */
