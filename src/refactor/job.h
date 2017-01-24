/*
 * job.h
 *
 *  Created on: Jan 24, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_REFACTOR_JOB_H_
#define SRC_REFACTOR_JOB_H_

#include <stdint.h>
#include <memory>

namespace is {

class Job;

std::shared_ptr<Job> createJob(int64_t seq_id, int64_t start, int64_t end, int64_t win_len);

class Job {
 public:
  ~Job() { };
  friend std::shared_ptr<Job> createJob(int64_t seq_id, int64_t start, int64_t end, int64_t win_len);

  int64_t end() const {
    return end_;
  }

  int64_t seq_id() const {
    return seq_id_;
  }

  int64_t start() const {
    return start_;
  }

  int64_t win_len() const {
    return win_len_;
  }

 private:
  explicit Job(int64_t ref_id, int64_t start, int64_t end, int64_t win_len);

  int64_t seq_id_;    // ID of the sequence from which the job is built.
  int64_t start_;     // 0-based start position for consensus
  int64_t end_;       // 0-based end position for consensus. Points to 1 base *after* the last inclusive base.
  int64_t win_len_;
};

} /* namespace is */

#endif /* SRC_REFACTOR_JOB_H_ */
