#ifndef INCLUDE_UNODE_ARENA
#define INCLUDE_UNODE_ARENA

#include <cstddef>
#include <vector>
#include <new>
#include <utility>
#include "unode.h"

// Block arena allocator for unode objects.
// Allocates unodes from contiguous blocks, freeing them all at once on
// destruction. Eliminates per-node malloc/free overhead in the hot
// utree copy path (branch-and-bound search).
//
// Safe because unode and all its members (neighbor_list, contracted_list,
// ints, bools) are trivially destructible — no per-object teardown needed.
class unode_arena {
	static constexpr std::size_t BLOCK_SIZE = 32;

	std::vector<void*> blocks_;
	std::size_t offset_ = BLOCK_SIZE;

public:
	unode_arena() = default;

	unode_arena(const unode_arena&) = delete;
	unode_arena& operator=(const unode_arena&) = delete;

	unode_arena(unode_arena&& other) noexcept
		: blocks_(std::move(other.blocks_)), offset_(other.offset_) {
		other.offset_ = BLOCK_SIZE;
	}

	unode_arena& operator=(unode_arena&& other) noexcept {
		if (this != &other) {
			destroy();
			blocks_ = std::move(other.blocks_);
			offset_ = other.offset_;
			other.offset_ = BLOCK_SIZE;
		}
		return *this;
	}

	~unode_arena() { destroy(); }

	template<typename... Args>
	unode* create(Args&&... args) {
		if (offset_ >= BLOCK_SIZE) {
			blocks_.push_back(
				::operator new(BLOCK_SIZE * sizeof(unode))
			);
			offset_ = 0;
		}
		void* slot = static_cast<char*>(blocks_.back()) +
			offset_ * sizeof(unode);
		++offset_;
		return ::new (slot) unode(std::forward<Args>(args)...);
	}

	friend void swap(unode_arena &first, unode_arena &second) noexcept {
		using std::swap;
		swap(first.blocks_, second.blocks_);
		swap(first.offset_, second.offset_);
	}

private:
	void destroy() {
		for (void* block : blocks_) {
			::operator delete(block);
		}
		blocks_.clear();
		offset_ = BLOCK_SIZE;
	}
};

#endif
