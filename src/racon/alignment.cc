/*
 * alignment.cc
 *
 *  Created on: Jan 26, 2017
 *      Author: isovic
 */

#include "alignment.h"
#include "utility/utility_general.h"

#include <thread>

namespace is {

int Alignment::AlignOverlap(const SingleSequence& query, const SingleSequence& target, const Overlap& overlap, int64_t overlap_id, int64_t win_size, int64_t win_ext, std::shared_ptr<SampledOverlap> sampled) {
//  std::thread::id thread_id = std::this_thread::get_id();
//  if (thread_id == 0) {
//    LOG_ALL("\rAligning overlap: %ld / %ld (%.2f\%), skipped %ld / %ld (%.2f\%)", i, overlaps.size(), 100.0f*((float) (i + 1)) / ((float) overlaps.size())
//            , num_skipped_overlaps, overlaps.size(), 100.0f*((float) (num_skipped_overlaps)) / ((float) overlaps.size()));
//  }

  std::shared_ptr<SingleSequence> seq(new SingleSequence());
  seq->CopyFrom(query);

  int64_t Astart = overlap.Astart();
  int64_t Aend = overlap.Aend();

  if (overlap.Brev()) {
    seq->ReverseComplement();
    Astart = overlap.Alen() - overlap.Aend();
    Aend = overlap.Alen() - overlap.Astart(); // There used to be a "- 1" here in the previous version of Racon and
                                              // the previous commits of the refactoring, but this is wrong. The alignment
                                              // of a sequence which should be perfect (the exact sequence used in the tiling path)
                                              // turned out to be by 1 longer than the query and target, having one 'I' at the
                                              // beginning and one 'D' at the end.
                                              // This demistyfies whether the END coordinate is non-inclusive or not (it is not inclusive!).
  }

  int64_t aln_start = 0, aln_end = 0, editdist = 0;
  std::vector<uint8_t> alignment;
  int rcaln = Alignment::AlignNW(seq->get_data() + Astart, (Aend - Astart),
                          target.get_data() + overlap.Bstart(), (overlap.Bend() - overlap.Bstart()),
                          &aln_start, &aln_end, &editdist, alignment);

  if (!rcaln) {

    // if (!overlap.Brev()) {
    //   printf ("Query:\n");
    //   for (int64_t i = Astart; i<Aend; i++) {
    //     printf ("%c", seq->get_data()[i]);
    //   }
    //   printf ("\n");
    //   printf ("\n");
    //   printf ("Target:\n");
    //   printf ("overlap.Bstart() = %ld\n", overlap.Bstart());
    //   printf ("overlap.Bend() = %ld\n", overlap.Bend());
    //   for (int64_t i = overlap.Bstart(); i<(overlap.Bend()); i++) {
    //     printf ("%c", target.get_data()[i]);
    //   }
    //   printf ("\n");
    //   printf ("\n");
    //   int64_t nm = 0, ni = 0, nd = 0, nx = 0;
    //   for (int64_t i=0; i<alignment.size(); i++) {
    //     printf ("%c", "=IDX"[alignment[i]]);
    //     nm += (alignment[i] == 0) ? 1 : 0;
    //     ni += (alignment[i] == 1) ? 1 : 0;
    //     nd += (alignment[i] == 2) ? 1 : 0;
    //     nx += (alignment[i] == 3) ? 1 : 0;
    //   }
    //   printf ("\n");
    //   printf ("nm = %ld, nx = %ld, ni = %ld, nd = %ld\n", nm, nx, ni, nd);
    //   printf ("e = %.2f\n", 100.0 * ((float) (nx + ni + nd)) / ((float) (nm + nx + ni + nd)));
    //   fflush(stdout);
    //   // exit(1);
    // }

    sampled->set(overlap, overlap_id, alignment, win_size, win_ext);

  } else {

    return 1;

  }

  return 0;
}

int SampleAlignment(const std::vector<unsigned char> &alignment, int64_t pos_on_ref) {
  return 0;
}

int Alignment::AlignNW(const int8_t *q, int64_t qlen, const int8_t *t,
                       int64_t tlen, int64_t* start, int64_t *end,
                       int64_t *eddist, std::vector<unsigned char> &alignment) {

  if (q == NULL || t == NULL || qlen <= 0 || tlen <= 0)
    return 1;

  int alphabet = 128;
  int score = 0;
  unsigned char* temp_aln = NULL;
  int aln_len = 0;

  int *positions = NULL;
  int num_positions = 0;
  int *start_locations = NULL;
  int found_k = 0;

  int myers_return_code = edlibCalcEditDistance((const unsigned char *) q, qlen,
                                                (const unsigned char *) t, tlen,
                                                alphabet, -1, EDLIB_MODE_NW,
                                                true, true, &score, &positions,
                                                &start_locations,
                                                &num_positions, &temp_aln,
                                                &aln_len, &found_k);

  if (myers_return_code == EDLIB_STATUS_ERROR ||
      num_positions == 0 || aln_len == 0) {
    return EDLIB_STATUS_ERROR;
  }

  int64_t ref_len = Alignment::CalcTargetLen(temp_aln, aln_len);

  *start = positions[0] - (ref_len - 1);
  *end = positions[0];
  *eddist = (int64_t) score;

  alignment.assign(temp_aln, (temp_aln + aln_len));

  if (positions)
    free(positions);

  if (start_locations)
    free(start_locations);

  if (temp_aln)
    free(temp_aln);

  return 0;
}

int32_t Alignment::CalcTargetLen(const unsigned char *aln, int32_t len) {
  if (aln == NULL || len == 0)
    return 0;

  int32_t tlen = 0;       // Target length.

  for (int32_t i = 0; i < len; i++) {
    if (aln[i] == EDLIB_EDOP_MATCH ||
        aln[i] == EDLIB_EDOP_MISMATCH ||
        aln[i] == EDLIB_EDOP_DELETE)
      tlen += 1;
  }

  return tlen;
}

} /* namespace is */
