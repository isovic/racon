/*
 * job.cc
 *
 *  Created on: Jan 24, 2017
 *      Author: Ivan Sovic
 */

#include "job.h"

namespace is {

std::shared_ptr<Job> createJob(int64_t seq_id, int64_t start, int64_t end, int64_t win_len) {
  return std::shared_ptr<Job>(new Job(seq_id, start, end, win_len));
}

Job::Job(int64_t seq_id, int64_t start, int64_t end, int64_t win_len)
         : seq_id_(seq_id), start_(start), end_(end), win_len_(win_len) {
}

}


