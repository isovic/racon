/*
 * mhap.cc
 *
 *  Created on: May 29, 2016
 *      Author: isovic
 */

#include <map>
#include <algorithm>
#include "libs/edlib.h"
#include "libs/edlibcigar.h"
#include "mhap.h"
#include <omp.h>
#include <vector>
#include <iostream>
#include <unordered_map>

int ParseOverlapLine(const std::string &line, OverlapFormat overlap_format, const std::map<std::string, int64_t> &qname_to_ids, const std::map<std::string, int64_t> &rname_to_ids, OldOverlapLine &overlap_line) {
  int rv = 0;
  if (overlap_format == kOverlapFormatPAF) {
    int rv = overlap_line.ParsePAF(line, qname_to_ids, rname_to_ids);
    if (rv) { return 1; }
  } else if (overlap_format == kOverlapFormatMHAP) {
    int rv = overlap_line.ParseMHAP(line);
    if (rv) { return 1; }
  } else {
    return 1;
  }
  return 0;
}

int TryAddingOverlap(const std::string &line, OverlapFormat overlap_format, const std::map<std::string, int64_t> &qname_to_ids, const std::map<std::string, int64_t> &rname_to_ids, float error_rate, std::unordered_map<int64_t, OldOverlapLine> &fmap) {
  // Parse the line and check if everything went fine.
  OldOverlapLine overlap_line;
  ParseOverlapLine(line, overlap_format, qname_to_ids, rname_to_ids, overlap_line);

  // Do simple filtering.
  if (overlap_line.CheckConstraints(error_rate)) {
    return 1;
  }

  // Look-up the A sequence in the map. Update if it's better, or add if
  // nothing yet exists for this A.
  auto it = fmap.find(overlap_line.Aid);
  if (it == fmap.end() || overlap_line.shared_minmers > it->second.shared_minmers) {
    fmap[overlap_line.Aid] = overlap_line;
  }
  return 0;
}

int ParseUniqueAndFilterErrors(const std::string &overlap_path, OverlapFormat overlap_format, const std::map<std::string, int64_t> &qname_to_ids, const std::map<std::string, int64_t> &rname_to_ids, float error_rate, std::vector<OldOverlapLine> &ret_overlaps) {
  std::unordered_map<int64_t, OldOverlapLine> fmap;     // Filtering map.

  if (overlap_path != "-") {
    std::ifstream infile(overlap_path);
    std::string line;
    while (std::getline(infile, line))
    {
      if (line.size() == 0) { continue; }
      std::istringstream iss(line);
      TryAddingOverlap(line, overlap_format, qname_to_ids, rname_to_ids, error_rate, fmap);
    }

    infile.close();

  } else {
    std::string line;
    while (std::getline(std::cin, line))
    {
      if (line.size() == 0) { break; }
      std::istringstream iss(line);
      TryAddingOverlap(line, overlap_format, qname_to_ids, rname_to_ids, error_rate, fmap);
    }
  }

  ret_overlaps.clear();
  for (auto it = fmap.begin(); it != fmap.end(); it++) {
    ret_overlaps.push_back(it->second);
  }

  return 0;
}

int ParseAndFilterErrors(const std::string &overlap_path, OverlapFormat overlap_format, const std::map<std::string, int64_t> &qname_to_ids, const std::map<std::string, int64_t> &rname_to_ids, float error_rate, std::vector<OldOverlapLine> &ret_overlaps) {
  ret_overlaps.clear();

  if (overlap_path != "-") {
    std::ifstream infile(overlap_path);
    std::string line;
    while (std::getline(infile, line))
    {
      if (line.size() == 0) { continue; }
      std::istringstream iss(line);

      // Parse the line and check if everything went fine.
      OldOverlapLine overlap_line;
      ParseOverlapLine(line, overlap_format, qname_to_ids, rname_to_ids, overlap_line);

      // Do simple filtering.
      if (!overlap_line.CheckConstraints(error_rate)) {
        ret_overlaps.push_back(overlap_line);
      }
    }

    infile.close();

  } else {
    std::string line;
    while (std::getline(std::cin, line))
    {
      if (line.size() == 0) { break; }
      std::istringstream iss(line);

      // Parse the line and check if everything went fine.
      OldOverlapLine overlap_line;
      ParseOverlapLine(line, overlap_format, qname_to_ids, rname_to_ids, overlap_line);

      // Do simple filtering.
      if (!overlap_line.CheckConstraints(error_rate)) {
        ret_overlaps.push_back(overlap_line);
      }
    }
  }

  return 0;
}



int FilterMHAP(const std::vector<OldOverlapLine> &overlaps_in, std::vector<OldOverlapLine> &overlaps_out, float error_rate) {
  std::map<int64_t, OldOverlapLine> fmap;     // Filtering map.

  for (int64_t i=0; i<overlaps_in.size(); i++) {
    if (!overlaps_in[i].CheckConstraints(error_rate)) {
      auto it = fmap.find(overlaps_in[i].Aid);
      if (it == fmap.end() || overlaps_in[i].shared_minmers > it->second.shared_minmers) {
        fmap[overlaps_in[i].Aid] = overlaps_in[i];
//      } else {
//        if (overlaps_in[i].shared_minmers > it->second.shared_minmers) {
//        }
      }
    }
  }
  overlaps_out.clear();
  for (auto it = fmap.begin(); it != fmap.end(); it++) {
    overlaps_out.push_back(it->second);
  }

//  std::sort(overlaps_out.begin(), overlaps_out.end(), [](const MHAPLine &a, const MHAPLine &b) { return a.Bstart < b.Bstart; });

  return 0;
}

int FilterMHAPErc(const std::vector<OldOverlapLine> &overlaps_in, std::vector<OldOverlapLine> &overlaps_out, float error_rate) {
  overlaps_out.clear();
  overlaps_out.reserve(overlaps_in.size());
  for (int64_t i=0; i<overlaps_in.size(); i++) {
    if (!overlaps_in[i].CheckConstraints(error_rate)) {
      overlaps_out.push_back(overlaps_in[i]);
    }
  }
  return 0;
}

//int DuplicateAndSwitch(const std::vector<OldOverlapLine> &overlaps_in, std::vector<OldOverlapLine> &overlaps_out) {
//  overlaps_out.clear();
//  overlaps_out.reserve(overlaps_in.size());
//  for (int64_t i=0; i<overlaps_in.size(); i++) {
//    overlaps_out.push_back(overlaps_in[i]);
//    OldOverlapLine o = overlaps_in[i];
//    o.Switch();
//    if (o.Arev) { o.ReverseComplement(); }
//    overlaps_out.push_back(o);
//  }
//  return 0;
//}

int DoReverseComplementing(std::vector<OldOverlapLine> &overlaps, SequenceFile &reads) {
  std::vector<bool> is_reversed(reads.get_sequences().size(), false);

  for (int64_t i=0; i<overlaps.size(); i++) {
    auto id = overlaps[i].Aid - 1;
    if (overlaps[i].Brev) {
      if (is_reversed[id] == false) {
        // Reverse the sequence.
        reads.get_sequences()[id]->ReverseComplement();
      }

      // Even though the sequence needs to be reversed only once, if there are more than one overlap,
      // reverse the strand in the coordinate space anyway.
      // Reverse the overlap itself.
      int64_t Astart = overlaps[i].Astart;
      overlaps[i].Astart = overlaps[i].Alen - overlaps[i].Aend;
      overlaps[i].Aend = overlaps[i].Alen - Astart - 1;
      overlaps[i].Brev = 0;

      // Mark reversed.
      is_reversed[id] = true;
    }
  }

  return 0;
}

int AlignOverlaps(const SequenceFile &refs, const SequenceFile &reads, const std::vector<OldOverlapLine> &overlaps, int32_t num_threads, SequenceFile &aligned, bool verbose_debug) {
  aligned.Clear();

  // Generate the SAM file header, for debugging.
  std::vector<std::string> sam_header;
  sam_header.push_back(FormatString("@HD\tVN:1.0\tSO:unknown"));
  for (int64_t i=0; i<refs.get_sequences().size(); i++) {
    const SingleSequence *ref = refs.get_sequences()[i];
    sam_header.push_back(FormatString("@SQ\tSN:%s\tLN:%ld", ref->get_header(), ref->get_data_length()));
  }
  aligned.set_file_header(sam_header);

//  // Generate the reverse complemented references, important for alignment.
//  SequenceFile rev_refs;
//  for (int64_t i=0; i<refs.get_sequences().size(); i++) {
//    const SingleSequence *ref = refs.get_sequences()[i];
//    SingleSequence *seq = new SingleSequence();
//
//    seq->InitDatafromAscii((int8_t *) ref->get_data(), ref->get_data_length());
//    seq->InitHeader((char *) ref->get_header(), ref->get_header_length());
//    seq->set_sequence_id(ref->get_sequence_id());
//    seq->set_sequence_absolute_id(ref->get_sequence_absolute_id());
//    seq->ReverseComplement();
//    rev_refs.AddSequence(seq, true);
//  }

//  fprintf (stderr, "\n");
  int64_t num_skipped_overlaps = 0;

  #pragma omp parallel for num_threads(num_threads) shared(aligned) schedule(dynamic, 1)
  for (int64_t i=0; i<overlaps.size(); i++) {
    int32_t thread_id = omp_get_thread_num();

    if (verbose_debug == true && thread_id == 0) {
      LOG_ALL("\rAligning overlap: %ld / %ld (%.2f\%), skipped %ld / %ld (%.2f\%)", i, overlaps.size(), 100.0f*((float) (i + 1)) / ((float) overlaps.size()), num_skipped_overlaps, overlaps.size(), 100.0f*((float) (num_skipped_overlaps)) / ((float) overlaps.size()));
      fflush(stderr);
    }

    auto &mhap = overlaps[i];
    OldOverlapLine omhap = mhap;

    // Get the read.
    const SingleSequence* read = reads.get_sequences()[omhap.Aid - 1];
    // Get the reference. It could be reversed.
    const SingleSequence* ref = NULL;
    ref = refs.get_sequences()[omhap.Bid - 1];

    std::string ref_name(ref->get_header(), ref->get_header_length());

    SingleSequence *seq = new SingleSequence();
    seq->InitAllFromAscii((char *) read->get_header(), read->get_header_length(),
                          (int8_t *) read->get_data(), (int8_t *) read->get_quality(), read->get_data_length(),
                          aligned.get_sequences().size(), aligned.get_sequences().size());
    if (omhap.Brev) {
      seq->ReverseComplement();
      omhap.Astart = mhap.Alen - mhap.Aend;
      omhap.Aend = mhap.Alen - mhap.Astart - 1;
    }

//    LOG_ALL(", len = %ld", seq->get_sequence_length());
//    if (omhap.Alen != seq->get_sequence_length()) {
//      printf ("\nWrong len! i = %ld, qname: '%s', len = %ld, omhap.Alen = %ld\n", i, seq->get_header(), seq->get_sequence_length(), omhap.Alen);
//      exit(1);
//    }

    int64_t aln_start = 0, aln_end = 0, editdist = 0;
    std::vector<unsigned char> alignment;
    int rcaln = EdlibNWWrapper(seq->get_data() + omhap.Astart, (omhap.Aend - omhap.Astart),
                            ref->get_data() + omhap.Bstart, (omhap.Bend - omhap.Bstart),
                            &aln_start, &aln_end, &editdist, alignment);

//    if (i == 1710 || i == 1712 || i == 1714 || i == 1722) {
//      printf ("\ni = %ld, qname: '%s', len = %ld\n", i, seq->get_header(), seq->get_sequence_length());
//      printf ("Read coords start: %ld, end: %ld, len: %ld:\n%s\n", omhap.Astart, omhap.Aend, omhap.Alen, GetSubstring((char *) (seq->get_data() + omhap.Astart), (omhap.Aend - omhap.Astart)).c_str());
//      printf ("Ref coords start: %ld, end: %ld, len: %ld:\n%s\n", omhap.Bstart, omhap.Bend, omhap.Blen, GetSubstring((char *) (ref->get_data() + omhap.Bstart), (omhap.Bend - omhap.Bstart)).c_str());
//      printf ("\n");
//      fflush(stdout);
//    }

    if (!rcaln) {
      std::string cigar_string;
      int rccig = edlibAlignmentToCigar(&alignment[0], alignment.size(), EDLIB_CIGAR_EXTENDED, cigar_string);

      SequenceAlignment aln;
      aln.SetCigarFromString(cigar_string);
      aln.set_pos(omhap.Bstart + aln_start + 1);
      aln.set_flag((int32_t) (16 * omhap.Brev));     // 0 if fwd and 16 if rev.
      aln.set_mapq(40);
      aln.set_rname(ref_name);

      // Remove insertions at front.
      for (int32_t j=0; j<aln.cigar().size(); j++) {
        if (aln.cigar()[j].count == 0) { continue; }
        if (aln.cigar()[j].op != 'D') { break; }
        aln.cigar()[j].count = 0;
      }
      for (int32_t j=0; j<aln.cigar().size(); j++) {
        if (aln.cigar()[j].count == 0) { continue; }
        if (aln.cigar()[j].op != 'I') { break; }
        aln.cigar()[j].op = 'S';
      }
      // Remove insertions at back.
      for (int32_t j=(aln.cigar().size()-1); j>=0; j--) {
        if (aln.cigar()[j].count == 0) { continue; }
        if (aln.cigar()[j].op != 'D') { break; }
        aln.cigar()[j].count = 0;
      }
      for (int32_t j=(aln.cigar().size()-1); j>=0; j--) {
        if (aln.cigar()[j].count == 0) { continue; }
        if (aln.cigar()[j].op != 'I') { break; }
        aln.cigar()[j].op = 'S';
      }

      // If the overlap does not cover the entire read (and most likely it does not).
      if (omhap.Astart > 0) {
        if (aln.cigar().size() > 0 && aln.cigar().front().op == 'S') { aln.cigar().front().count += omhap.Astart; }
        else { CigarOp new_op; new_op.op = 'S'; new_op.count = omhap.Astart; aln.cigar().insert(aln.cigar().begin(), new_op); }
      }
      if ((omhap.Aend) < (omhap.Alen)) {
        if (aln.cigar().size() > 0 && aln.cigar().back().op == 'S') { aln.cigar().back().count += (omhap.Alen - omhap.Aend); }
        else { CigarOp new_op; new_op.op = 'S'; new_op.count = (omhap.Alen - omhap.Aend); aln.cigar().insert(aln.cigar().end(), new_op); }
      }

      aln.RecalcCigarPositions();

      int64_t m=0, x=0, eq=0, ins=0, del=0;
      for (int32_t j=0; j<aln.cigar().size(); j++) {
        if (aln.cigar()[j].op == 'M') { m += aln.cigar()[j].count; }
        if (aln.cigar()[j].op == '=') { eq += aln.cigar()[j].count; }
        if (aln.cigar()[j].op == 'X') { x += aln.cigar()[j].count; }
        if (aln.cigar()[j].op == 'I') { ins += aln.cigar()[j].count; }
        if (aln.cigar()[j].op == 'D') { del += aln.cigar()[j].count; }
      }
      aln.optional().push_back(FormatString("X1:Z:equal=%ld_x=%ld_ins=%ld_del=%ld", eq, x, ins, del));

      seq->InitAlignment(aln);

      #pragma omp critical
      aligned.AddSequence(seq, true);

    } else {
      if (seq) {
        delete seq;
        seq = NULL;
      }
      num_skipped_overlaps += 1;
    }

  }

  if (verbose_debug == true) {
    fprintf (stderr, "\n");
    fflush(stderr);
  }

//  FILE *fp = fopen("temp/debug.sam", "w");
//  for (int32_t i=0; i<aligned.get_file_header().size(); i++) {
//    fprintf (fp, "%s\n", aligned.get_file_header()[i].c_str());
//  }
//  for (int64_t i=0; i<aligned.get_sequences().size(); i++) {
//    std::string line = aligned.get_sequences()[i]->MakeSAMLine();
//    fprintf (fp, "%s\n", line.c_str());
//  }
//  fclose(fp);

  return 0;
}

int64_t CalculateReconstructedLength(unsigned char *alignment, int alignmentLength) {
  if (alignment == NULL || alignmentLength == 0)
      return 0;

  int64_t length = 0;
  std::stringstream ss;

  for (int i=0; i<alignmentLength; i++) {
    if (alignment[i] == EDLIB_EDOP_MATCH || alignment[i] == EDLIB_EDOP_MISMATCH || alignment[i] == EDLIB_EDOP_DELETE)
      length += 1;
  }

  return length;
}

int EdlibNWWrapper(const int8_t *read_data, int64_t read_length,
                   const int8_t *reference_data, int64_t reference_length,
                   int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                   int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment) {

  if (read_data == NULL || reference_data == NULL || read_length <= 0 || reference_length <= 0)
    return 1;

  int alphabet_length = 128;
  int score = 0;
  unsigned char* alignment = NULL;
  int alignment_length = 0;

  int *positions = NULL;
  int num_positions = 0;
  int *start_locations = NULL;
  int found_k = 0;

  int myers_return_code = edlibCalcEditDistance((const unsigned char *) read_data, read_length,
                        (const unsigned char *) reference_data, reference_length,
                        alphabet_length, -1, EDLIB_MODE_NW, true, true, &score, &positions, &start_locations, &num_positions,
                        &alignment, &alignment_length, &found_k);

  if (myers_return_code == EDLIB_STATUS_ERROR || num_positions == 0 || alignment_length == 0) {
    return EDLIB_STATUS_ERROR;
  }

  int64_t reconstructed_length = CalculateReconstructedLength(alignment, alignment_length);

  *ret_alignment_position_start = positions[0] - (reconstructed_length - 1);
//  *ret_alignment_position_start = start_locations[0];
  *ret_alignment_position_end = positions[0];
  *ret_edit_distance = (int64_t) score;
  ret_alignment.assign(alignment, (alignment + alignment_length));

  if (positions)
    free(positions);

  if (start_locations)
    free(start_locations);

  if (alignment)
      free(alignment);

  return 0;
}
