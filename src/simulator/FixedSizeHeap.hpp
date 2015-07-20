#ifndef FIXEDSIZEHEAP_HPP
#define FIXEDSIZEHEAP_HPP

#include <stdexcept>
#include <vector>

// The heap size is fixed on construction, designed to conserve memory.
// The heap elements are stored in an array.
// The children of node k are nodes 2k+1 and 2k+2.
// The heap property requires that the tree have a specific shape:
// all the levels are full except for possibly the lowest level,
// and the lowest level is filled in from the left.
// Also, each node have a value no larger than its children.
// This implementation based on Chapter 12 of Jon Bentley's "Programming Pearls", 1989.

// Attempting to insert an element into a full array causes an exception.
// The user is responsible for determining whether the array is full
// before attempting an insert.

// By John D. Cook, http://www.johndcook.com

template <typename TStorageType = float>
class FixedSizeHeap
{
public:

    FixedSizeHeap(int size = 0)
    {
        Initialize(size);
    }

    void Initialize(int size)
    {
        m_heap.resize(size);
        m_lastIndex = -1;
    }

    // Return the effective size, the number of slots in the heap that are used and not just allocated.
    int EffectiveSize()
    {
        return m_lastIndex + 1;
    }

    // Look at the smallest element but leave it alone
    TStorageType InspectMinElement()
    {
        return m_heap[0];
    }

    // Replace the smallest element with the input argument
    void ReplaceMinElement(TStorageType x)
    {
        m_heap[0] = x;
        SiftDown();
    }

    // Insert an element into the heap maintaining the heap property
    void Insert(TStorageType x)
    {
        m_lastIndex++;
        if (m_lastIndex >= static_cast<int>(m_heap.size()))
        {
            m_lastIndex--;
			throw std::length_error("Heap is full");
        }
        m_heap[m_lastIndex] = x;
        SiftUp();
    }

protected:
	std::vector<TStorageType> m_heap;
    int m_lastIndex;

	// The tester class is a friend so it can test protected methods of this class.
    friend class FixedSizeHeapTester;

    // Move the effective last element of the array into the
    // necessary position to restore the heap property.
    // This moves the element "up" the tree, envisioning the root at the top.
    void SiftUp()
    {
        for (int i = m_lastIndex; i > 0; /* no decrement */)
        {
            // Find the parent node.
            int parent = (i - 1)/2;  // deliberate integer division

            // If the heap property holds, then stop.
            if (m_heap[parent]<= m_heap[i])
                break;

            // Otherwise, swap up the tree until it is.
            TStorageType temp = m_heap[parent]; m_heap[parent] = m_heap[i]; m_heap[i] = temp;
            i = parent;
       }
    }

    // Move the first element of the array into the
    // necessary position to restore the heap property.
    // This moves the element "down" the tree, envisioning the root at the top.
    void SiftDown()
    {
        int allocatedSize = static_cast<int>(m_heap.size());
        for (int i = 0; ; )
        {
            int child = 2*i + 1;  // "child" is the left child of node i
            if (child >= allocatedSize) // node i has no left child within allocated space
                break;
            int rightChild = child + 1;
            if (rightChild < allocatedSize) // there exists a right child within allocated space
            {
                // Make "child" the least child of i.
                if (m_heap[rightChild] < m_heap[child])
                    child = rightChild;
			}

            // If heap property holds, then stop.
            if (m_heap[i] <= m_heap[child])
                break;

            // Otherwise, swap down the tree.
            TStorageType temp = m_heap[i]; m_heap[i] = m_heap[child]; m_heap[child] = temp;
            i = child;
        }
    }
};

#endif
