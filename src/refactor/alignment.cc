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
    Aend = overlap.Alen() - overlap.Astart() - 1;
  }

  int64_t aln_start = 0, aln_end = 0, editdist = 0;
  std::vector<uint8_t> alignment;
  int rcaln = Alignment::AlignNW(seq->get_data() + Astart, (Aend - Astart),
                          target.get_data() + overlap.Bstart(), (overlap.Bend() - overlap.Bstart()),
                          &aln_start, &aln_end, &editdist, alignment);

  if (!rcaln) {
    sampled->set(overlap, overlap_id, alignment, win_size, win_ext);

  } else {
    return 1;

  }

  return 0;
}

//int Alignment::AlignOverlaps(const SequenceFile &refs, const SequenceFile &reads, const std::vector<OldOverlapLine> &overlaps, int32_t num_threads, SequenceFile &aligned, bool verbose_debug) {
//  aligned.Clear();
//
//  // Generate the SAM file header, for debugging.
//  std::vector<std::string> sam_header;
//  sam_header.push_back(FormatString("@HD\tVN:1.0\tSO:unknown"));
//  for (int64_t i=0; i<refs.get_sequences().size(); i++) {
//    const SingleSequence *ref = refs.get_sequences()[i];
//    sam_header.push_back(FormatString("@SQ\tSN:%s\tLN:%ld", ref->get_header(), ref->get_data_length()));
//  }
//  aligned.set_file_header(sam_header);
//
//  int64_t num_skipped_overlaps = 0;
//
//  #pragma omp parallel for num_threads(num_threads) shared(aligned) schedule(dynamic, 1)
//  for (int64_t i=0; i<overlaps.size(); i++) {
//    int32_t thread_id = omp_get_thread_num();
//
//    if (verbose_debug == true && thread_id == 0) {
//      LOG_ALL("\rAligning overlap: %ld / %ld (%.2f\%), skipped %ld / %ld (%.2f\%)", i, overlaps.size(), 100.0f*((float) (i + 1)) / ((float) overlaps.size()), num_skipped_overlaps, overlaps.size(), 100.0f*((float) (num_skipped_overlaps)) / ((float) overlaps.size()));
//      fflush(stderr);
//    }
//
//    auto &mhap = overlaps[i];
//    OldOverlapLine omhap = mhap;
//
//    // Get the read.
//    const SingleSequence* read = reads.get_sequences()[omhap.Aid - 1];
//    // Get the reference. It could be reversed.
//    const SingleSequence* ref = NULL;
//    ref = refs.get_sequences()[omhap.Bid - 1];
//
//    std::string ref_name(ref->get_header(), ref->get_header_length());
//
//    SingleSequence *seq = new SingleSequence();
//    seq->InitAllFromAscii((char *) read->get_header(), read->get_header_length(),
//                          (int8_t *) read->get_data(), (int8_t *) read->get_quality(), read->get_data_length(),
//                          aligned.get_sequences().size(), aligned.get_sequences().size());
//    if (omhap.Brev) {
//      seq->ReverseComplement();
//      omhap.Astart = mhap.Alen - mhap.Aend;
//      omhap.Aend = mhap.Alen - mhap.Astart - 1;
//    }
//
//    int64_t aln_start = 0, aln_end = 0, editdist = 0;
//    std::vector<unsigned char> alignment;
//    int rcaln = EdlibNWWrapper(seq->get_data() + omhap.Astart, (omhap.Aend - omhap.Astart),
//                            ref->get_data() + omhap.Bstart, (omhap.Bend - omhap.Bstart),
//                            &aln_start, &aln_end, &editdist, alignment);
//
//    if (!rcaln) {
//      std::string cigar_string;
//      int rccig = edlibAlignmentToCigar(&alignment[0], alignment.size(), EDLIB_CIGAR_EXTENDED, cigar_string);
//
//      SequenceAlignment aln;
//      aln.SetCigarFromString(cigar_string);
//      aln.set_pos(omhap.Bstart + aln_start + 1);
//      aln.set_flag((int32_t) (16 * omhap.Brev));     // 0 if fwd and 16 if rev.
//      aln.set_mapq(40);
//      aln.set_rname(ref_name);
//
//      // Remove insertions at front.
//      for (int32_t j=0; j<aln.cigar().size(); j++) {
//        if (aln.cigar()[j].count == 0) { continue; }
//        if (aln.cigar()[j].op != 'D') { break; }
//        aln.cigar()[j].count = 0;
//      }
//      for (int32_t j=0; j<aln.cigar().size(); j++) {
//        if (aln.cigar()[j].count == 0) { continue; }
//        if (aln.cigar()[j].op != 'I') { break; }
//        aln.cigar()[j].op = 'S';
//      }
//      // Remove insertions at back.
//      for (int32_t j=(aln.cigar().size()-1); j>=0; j--) {
//        if (aln.cigar()[j].count == 0) { continue; }
//        if (aln.cigar()[j].op != 'D') { break; }
//        aln.cigar()[j].count = 0;
//      }
//      for (int32_t j=(aln.cigar().size()-1); j>=0; j--) {
//        if (aln.cigar()[j].count == 0) { continue; }
//        if (aln.cigar()[j].op != 'I') { break; }
//        aln.cigar()[j].op = 'S';
//      }
//
//      // If the overlap does not cover the entire read (and most likely it does not).
//      if (omhap.Astart > 0) {
//        if (aln.cigar().size() > 0 && aln.cigar().front().op == 'S') { aln.cigar().front().count += omhap.Astart; }
//        else { CigarOp new_op; new_op.op = 'S'; new_op.count = omhap.Astart; aln.cigar().insert(aln.cigar().begin(), new_op); }
//      }
//      if ((omhap.Aend) < (omhap.Alen)) {
//        if (aln.cigar().size() > 0 && aln.cigar().back().op == 'S') { aln.cigar().back().count += (omhap.Alen - omhap.Aend); }
//        else { CigarOp new_op; new_op.op = 'S'; new_op.count = (omhap.Alen - omhap.Aend); aln.cigar().insert(aln.cigar().end(), new_op); }
//      }
//
//      aln.RecalcCigarPositions();
//
//      int64_t m=0, x=0, eq=0, ins=0, del=0;
//      for (int32_t j=0; j<aln.cigar().size(); j++) {
//        if (aln.cigar()[j].op == 'M') { m += aln.cigar()[j].count; }
//        if (aln.cigar()[j].op == '=') { eq += aln.cigar()[j].count; }
//        if (aln.cigar()[j].op == 'X') { x += aln.cigar()[j].count; }
//        if (aln.cigar()[j].op == 'I') { ins += aln.cigar()[j].count; }
//        if (aln.cigar()[j].op == 'D') { del += aln.cigar()[j].count; }
//      }
//      aln.optional().push_back(FormatString("X1:Z:equal=%ld_x=%ld_ins=%ld_del=%ld", eq, x, ins, del));
//
//      seq->InitAlignment(aln);
//
//      #pragma omp critical
//      aligned.AddSequence(seq, true);
//
//    } else {
//      if (seq) {
//        delete seq;
//        seq = NULL;
//      }
//      num_skipped_overlaps += 1;
//    }
//
//  }
//
//  if (verbose_debug == true) {
//    fprintf (stderr, "\n");
//    fflush(stderr);
//  }
//
////  FILE *fp = fopen("temp/debug.sam", "w");
////  for (int32_t i=0; i<aligned.get_file_header().size(); i++) {
////    fprintf (fp, "%s\n", aligned.get_file_header()[i].c_str());
////  }
////  for (int64_t i=0; i<aligned.get_sequences().size(); i++) {
////    std::string line = aligned.get_sequences()[i]->MakeSAMLine();
////    fprintf (fp, "%s\n", line.c_str());
////  }
////  fclose(fp);
//
//  return 0;
//}

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
