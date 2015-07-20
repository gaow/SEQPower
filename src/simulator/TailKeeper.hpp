#ifndef TAILKEEPER_HPP
#define TAILKEEPER_HPP

#include "FixedSizeHeap.hpp"

// This class is used for finding percentiles without storing an entire set of data.
// As samples enter, only the tails are preserved: if we know a value cannot be in the
// left or right tail, it is ignored.  Thus the storage size is determined by the tail
// sizes rather than the number of samples.  For example, the left and right 5% tails
// can be stored with only 10% of the storage necessary to save and sort all samples.

// The input and storage types may differ.  The motivating application would enter
// doubles but saving them as floats to conserve memory.

// John D. Cook, http://www.johndcook.com


template <typename TInputType = double, typename TStorageType = float>
class TailKeeper
{
public:
    TailKeeper(int leftTailSize = 0, int rightTailSize = 0)
    {
        Initialize(leftTailSize, rightTailSize);
    }

    void Initialize(int leftTailSize, int rightTailSize)
    {
        m_leftTailSize  = leftTailSize;
        m_rightTailSize = rightTailSize;

        m_leftTail.Initialize(leftTailSize);
        m_rightTail.Initialize(rightTailSize);
    }

    void AddSample(TInputType x)
    {
        TStorageType castedx = static_cast<TStorageType>(x);

        // Save largest elements, i.e. right tail
        if (m_rightTail.EffectiveSize() < m_rightTailSize)
            m_rightTail.Insert(castedx);
        else if (m_rightTail.InspectMinElement() < castedx)
            m_rightTail.ReplaceMinElement(castedx);

        // Save smallest elements, i.e. left tail.
        // NB: This one is slightly trickier.
        // The heap is designed to designed to maintain the heap property in terms of
        // MINIMUMS, but we want to use it to keep track of the MAXIMUMS of the left tail.
        // Therefore we insert the NEGATIVES of the elements.
        castedx = -castedx;
        if (m_leftTail.EffectiveSize() < m_leftTailSize)
            m_leftTail.Insert(castedx);
        else if (m_leftTail.InspectMinElement() < castedx)
            m_leftTail.ReplaceMinElement(castedx);
    }

    TInputType GetMaxLeftTail()
    {
        // See explanation in AddSample for why there is a negative sign here.
        return static_cast<TInputType>( -m_leftTail.InspectMinElement() );
    }

    TInputType GetMinRightTail()
    {
       return static_cast<TInputType>( m_rightTail.InspectMinElement() );
    }


protected:
    FixedSizeHeap<TStorageType> m_leftTail;
    FixedSizeHeap<TStorageType> m_rightTail;
    int m_leftTailSize;
    int m_rightTailSize;

};

#endif
