//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#ifndef KRATOS_KEY_HASH_H_INCLUDED
#define KRATOS_KEY_HASH_H_INCLUDED

namespace Kratos 
{
///@addtogroup KratosCore
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the hash type
    typedef std::size_t HashType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

    /**
     * @brief This method creates an "unique" hash for the input value
     * @tparam TClassType The type of class to be hashed
     * @param Seed This is the seed used to create the hash
     * @param Value This is the value to be hashed
     */
    template <class TClassType>
    inline void HashCombine(
        HashType& Seed,
        const TClassType& Value
        )
    {
        std::hash<TClassType> hasher;
        Seed ^= hasher(Value) + 0x9e3779b9 + (Seed<<6) + (Seed>>2);
    }
    
    /**
     * @brief This method combines hash until it reaches the last class in order to obtain a corresponding seed
     * @tparam TClassType The type of class to be hashed
     * @param First The first class to be compared
     * @param Last The last class to be compared
     * @return The resulting seed
     */
    template <class TClassType>
    inline HashType HashRange(
        TClassType First, 
        TClassType Last
        )
    {
        HashType seed = 0;

        while (First!=Last) {
            HashCombine(seed, *First);
            ++First;
        }
        
        return seed;
    }
    
    /**
     * @brief This is a key comparer of general pourpose between two classes
     * @tparam TClassType The type of class to be hashed
     */
    template<class TClassType>
    struct KeyComparorRange
    {
        /**
         * @brief This is the () operator
         * @param first The first class to be compared
         * @param second The second class to be compared
         */
        bool operator()(
            const TClassType& first, 
            const TClassType& second
            ) const
        {
            if(first.size() != second.size()) {
                return false;
            }

            auto it_first = first.begin();
            auto it_second = second.begin();

            while(it_first != first.end()) { // NOTE: We already checked that are same size
                if(*it_first != *it_second) {
                    return false;
                }
                if(it_first != first.end()) {
                    ++it_first;
                    ++it_second;
                }
            }

            return true;
        }
    };
    
    /**
     * @brief This is a hasher of general pourpose
     * @tparam TClassType The type of class to be hashed
     */
    template<class TClassType>
    struct KeyHasherRange
    {
        /**
         * @brief This is the () operator
         * @param rRange The shared pointer to be hashed
         * @return The corresponding hash
         */
        HashType operator()(const TClassType& rRange) const
        {
            return HashRange(rRange.begin(), rRange.end());
        }
    };
    
    /**
     * @brief This is a hasher for shared pointers
     * @tparam TSharedPointer The type of shared pointer to be hashed
     */
    template<class TSharedPointer>
    struct SharedPointerHasher
    {
        /**
         * @brief This is the () operator
         * @param pPointer The shared pointer to be hashed
         * @return The corresponding hash
         */
        HashType operator()(const TSharedPointer& pPointer) const
        {
            return reinterpret_cast<HashType>(pPointer.get());
        }
    };
    
    /**
     * @brief This is a key comparer between two shared pointers
     * @tparam TSharedPointer The type of shared pointer to be compared
     */
    template<class TSharedPointer>
    struct SharedPointerComparator
    {
        /**
         * @brief This is the () operator
         * @param first The first class to be compared
         * @param second The second class to be compared
         */
        bool operator()(
            const TSharedPointer& first, 
            const TSharedPointer& second
            ) const
        {
            return first.get() == second.get();
        }
    };

    /**
     * @brief This is a hasher between two vectors of indexes
     * @tparam TVectorIndex The type of vector indexes to be compared
     */
    template<class TVectorIndex>
    struct VectorIndexHasher
    {
        /**
         * @brief This is the () operator
         * @param k The vector of indexes to be hashed
         * @return The corresponding hash
         */
        HashType operator()(const TVectorIndex& k) const
        {
            return HashRange(k.begin(), k.end());
        }
    };

    /**
     * @brief This is a key comparer between two vectors of indexes
     * @tparam TVectorIndex The type of vector indexes to be compared
     */
    template<class TVectorIndex>
    struct VectorIndexComparor
    {
        /**
         * @brief This is the () operator
         * @param lhs The first class to be compared
         * @param rhs The second class to be compared
         */
        bool operator()(const TVectorIndex& lhs, const TVectorIndex& rhs) const
        {
            if(lhs.size() != rhs.size())
                return false;

            for(IndexType i = 0; i < lhs.size(); i++) {
                if(lhs[i] != rhs[i]) return false;
            }

            return true;
        }
    };

    
///@}
///@name Kratos Classes
///@{

} // namespace Kratos.
#endif
