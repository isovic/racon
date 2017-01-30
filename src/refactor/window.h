/*
 * window.h
 *
 *  Created on: Jan 29, 2017
 *      Author: isovic
 */

#ifndef SRC_REFACTOR_WINDOW_H_
#define SRC_REFACTOR_WINDOW_H_

#include "types.h"

namespace is {

/** A container object holding info about a single query overlapping
 *  a particular window.
*/
class WindowEntry {
 public:
  WindowEntry() : overlap_id_(-1) { }
  WindowEntry(int64_t overlap_id, int64_t qstart, int64_t qend, int64_t tstart, int64_t tend) :
   			overlap_id_(overlap_id), query_.start(qstart), query_.end(qend), target_.start(tstart), target_.end(tend) { }
  ~WindowEntry() { }

  int64_t overlap_id() const {
  	return overlap_id_;
  }
  Range query() const {
  	return query;
  }
  Range target() const {
  	return target;
  }

 private:
  int64_t overlap_id_ = - 1;
  Range query_, target_;					// Query and target start and end locations.
};

class Window {
 public:
  Window();
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

 private:
//  std::vector<
  std::vector<WindowEntry> entries_;
};

} /* namespace is */

#endif /* SRC_REFACTOR_WINDOW_H_ */
