/*!
* @file cudaaligner.hpp
 *
 * @brief CUDA aligner class header file
 */
#include "cudaaligner/cudaaligner.hpp"
#include "cudaaligner/aligner.hpp"
#include "cudaaligner/alignment.hpp"

#include "overlap.hpp"
#include "sequence.hpp"

#include <vector>
#include <atomic>

namespace racon {

class CUDABatchAligner;
std::unique_ptr<CUDABatchAligner> createCUDABatchAligner(uint32_t max_query_size, uint32_t max_target_size, uint32_t max_alignments, uint32_t device_id);

class CUDABatchAligner
{
    public:
        virtual ~CUDABatchAligner();

        /**
         * @brief Add a new overlap to the batch.
         *
         * @param[in] window   : The overlap to add to the batch.
         * @param[in] sequences: Reference to a database of sequences.
         *
         * @return True if overlap could be added to the batch.
         */
        virtual bool addOverlap(Overlap* overlap, std::vector<std::unique_ptr<Sequence>>& sequences);

        /**
         * @brief Checks if batch has any overlaps to process.
         *
         * @return Trie if there are overlaps in the batch.
         */
        virtual bool hasOverlaps() const {
            return overlaps_.size() > 0;
        };

        /**
         * @brief Runs batched alignment of overlaps on GPU.
         *
         */
        virtual void alignAll();

        /**
         * @brief Find breaking points in alignments.
         *
         */
        virtual void find_breaking_points(uint32_t window_length);

        /**
         * @brief Resets the state of the object, which includes
         *        resetting buffer states and counters.
         */
        virtual void reset();

        /**
         * @brief Get batch ID.
         */
        uint32_t getBatchID() const { return bid_; }

        // Builder function to create a new CUDABatchAligner object.
        friend std::unique_ptr<CUDABatchAligner>
        createCUDABatchAligner(uint32_t max_query_size, uint32_t max_target_size, uint32_t max_alignments, uint32_t device_id);

    protected:
        CUDABatchAligner(uint32_t max_query_size, uint32_t max_target_size, uint32_t max_alignments, uint32_t device_id);
        CUDABatchAligner(const CUDABatchAligner&) = delete;
        const CUDABatchAligner& operator=(const CUDABatchAligner&) = delete;

        std::unique_ptr<genomeworks::cudaaligner::Aligner> aligner_;

        std::vector<Overlap*> overlaps_;

        // Static batch count used to generate batch IDs.
        static std::atomic<uint32_t> batches;

        // Batch ID.
        uint32_t bid_ = 0;

        // CUDA stream for batch.
        cudaStream_t stream_;
};

}
