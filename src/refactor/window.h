/*
 * window.h
 *
 *  Created on: Jan 29, 2017
 *      Author: isovic
 */

#ifndef SRC_REFACTOR_WINDOW_H_
#define SRC_REFACTOR_WINDOW_H_

namespace is {

struct WindowEntry {
  int64_t query_id;

};

class Window {
 public:
  Window();
  ~Window();

  /** Adds info for a part of a sequence overlapping this particular window.
   * @overlap_id ID of the overlap used to generate this info.
   * @qstart Query
   */
  void add(int64_t overlap_id, int64_t qstart, int64_t qend, int64_t rstart, int64_t rend);

 private:
//  std::vector<
  std::vector<WindowEntry> entries_;
};

} /* namespace is */

#endif /* SRC_REFACTOR_WINDOW_H_ */
