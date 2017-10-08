/*
 * window.cc
 *
 *  Created on: Jan 29, 2017
 *      Author: isovic
 */

#include "window.h"

namespace is {

Window::Window(int64_t target_id) : target_id_(target_id), start_(0), end_(0) {

}

Window::~Window() {

}

void Window::SortEntries(bool skip_first) {
  if (entries_.size() == 0) { return; }

  auto it = entries_.begin();
  if (skip_first) {
    it++;
  }

  std::sort(it, entries_.end(), [](const WindowEntry& a, const WindowEntry& b) { return a.target().len() > b.target().len(); });

}



} /* namespace is */
