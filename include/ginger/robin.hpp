#pragma once

#include <vector>

namespace fun {

    namespace detail {
        /**
         * @brief A node in the Robin Hood hashing linked list
         *
         * @tparam T The type of the key stored in the node
         */
        template <typename T> struct RobinSlNode {
            RobinSlNode *next;
            T key;
        };

        /**
         * @brief Iterator for the Robin Hood linked list
         *
         * @tparam T The type of elements being iterated
         */
        template <typename T> struct RobinIterator {
            const RobinSlNode<T> *cur;
            auto operator!=(const RobinIterator &other) const -> bool { return cur != other.cur; }
            auto operator==(const RobinIterator &other) const -> bool { return cur == other.cur; }
            auto operator++() -> RobinIterator & {
                cur = cur->next;
                return *this;
            }
            auto operator*() const -> const T & { return cur->key; }
        };

        /**
         * @brief Wrapper to make the linked list iterable
         *
         * @tparam T The type of elements in the list
         */
        template <typename T> struct RobinIterableWrapper {
            const detail::RobinSlNode<T> *node;
            // const Robin<T> *rr;
            // T from_part;
            auto begin() const -> RobinIterator<T> { return RobinIterator<T>{node->next}; }
            auto end() const -> RobinIterator<T> { return RobinIterator<T>{node}; }
            // auto size() const -> size_t { return rr->cycle.size() - 1; }
        };
    }  // namespace detail

    /**
     * @brief Robin Hood iterator for cycle iteration
     *
     * A data structure that creates a circular linked list allowing iteration
     * over all elements except a specified one (the "exclude" operation).
     *
     * @tparam T The type of indices in the cycle
     */
    template <typename T> struct Robin {
        std::vector<detail::RobinSlNode<T>> cycle;

        /**
         * @brief Construct a new Robin cycle
         *
         * Creates a circular linked list containing all indices from 0 to num_parts-1.
         *
         * @param num_parts The number of elements in the cycle
         */
        explicit Robin(T num_parts) : cycle(num_parts) {
            auto *slptr = &this->cycle[num_parts - 1];
            auto k = T(0);
            for (auto &sl : this->cycle) {
                sl.key = k;
                slptr->next = &sl;
                slptr = slptr->next;
                ++k;
            }
        }

        /**
         * @brief Get an iterable that excludes a specific element
         *
         * Returns an iterable that allows iterating over all elements in the cycle
         * except the one at the specified index.
         *
         * @param from_part The index to exclude from iteration
         * @return detail::RobinIterableWrapper<T> An iterable over the remaining elements
         */
        auto exclude(T from_part) const -> detail::RobinIterableWrapper<T> {
            return detail::RobinIterableWrapper<T>{&this->cycle[from_part]};
        }
    };

}  // namespace fun
