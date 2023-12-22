#include <iostream>
#include <vector>
#include <algorithm>

/*
 * Insertion sort -- O(N^2)
 * @param elements: a reference to vector of elements to be sorted
 *
 * Done in-place, so we return nothing
 */
template<typename T>
void insertion_sort(std::vector<T> &elements) {
    if (elements.empty() || elements.size() <= 1) return;

    for (size_t i = 1; i < elements.size(); i++) {
        size_t j = i;
        while (j > 0 && elements[j - 1] > elements[j]) {
            std::swap(elements[j - 1], elements[j]);
            j--;
        }
    }
}

/*
 * Bubble sort -- O(N^2)
 * @param elements: a reference to vector of elements to be sorted
 *
 * Done in-place, so we return nothing
 */
template<typename T>
void bubble_sort(std::vector<T> &elements) {
    if (elements.empty() || elements.size() == 1) return;

    for (size_t i = 0; i < elements.size(); i++) {
        for (size_t j = 0; j < elements.size() - 1 - i; j++) {
            if (elements[j] > elements[j + 1]) {
                std::swap(elements[j], elements[j + 1]);
            }
        }
    }
}

/*
 * Selection sort -- O(N^2)
 * @param elements: a reference to vector of elements to be sorted
 *
 * Done in-place, so we return nothing
 */
template<typename T>
void selection_sort(std::vector<T> &elements) {
    if (elements.empty() || elements.size() == 1) return;

    for (size_t i = 0; i < elements.size(); i++) {
        size_t min_index_position = i;

        for (size_t j = i; j < elements.size(); j++) {
            if (elements[j] < elements[min_index_position]) {
                min_index_position = j;
            }
        }

        std::swap(elements[i], elements[min_index_position]);
    }
}

/*
 * Quick sort -- O(N log N) average case, worst case O(N^2)
 * @param elements: a reference to vector of elements to be sorted
 *
 * Not in place, so we can understand the algorithm better
 */
template<typename T>
std::vector<T> naive_quick_sort(const std::vector<T> &elements) { // NOLINT
    if (elements.empty() || elements.size() == 1) return elements;

    auto pivot =  elements[0];
    auto vec_left = std::vector<T>();
    auto vec_right = std::vector<T>();

    for (size_t i = 1; i < elements.size(); i++) {
        if (elements[i] < pivot) {
            vec_left.emplace_back(elements[i]);
        } else {
            vec_right.emplace_back(elements[i]);
        }
    }

    auto sorted_left = naive_quick_sort(vec_left);
    auto sorted_right = naive_quick_sort(vec_right);

    std::vector<T> result;
    result.reserve(sorted_left.size() + 1 + sorted_right.size());

    // the equivalent of [sorted_left] + pivot + [sorted_right]
    result.insert(result.end(), sorted_left.begin(), sorted_left.end());
    result.push_back(pivot);
    result.insert(result.end(), sorted_right.begin(), sorted_right.end());

    return result;
}

/*
 * Quick sort -- O(N log N) average case, worst case O(N^2)
 * @param elements: a reference to vector of elements to be sorted
 *
 * Done in-place, so we return nothing
 */
template <typename T>
void quick_sort_partition(std::vector<T>& elements, int start, int end) { // NOLINT
    if (start >= end) return;

    T pivot = elements[start];
    int i = start + 1;
    int j = end;

    while (i <= j) {
        while (i <= end && elements[i] <= pivot) i++;
        while (j > start && elements[j] > pivot) j--;

        if (i < j) {
            std::swap(elements[i], elements[j]);
        }
    }

    // Move the pivot element to its final position
    int pivot_index = j;
    std::swap(elements[start], elements[pivot_index]);

    quick_sort_partition(elements, start, pivot_index - 1);
    quick_sort_partition(elements, pivot_index + 1, end);
}


template <typename T>
void quick_sort(std::vector<T> &elements) {
    if (elements.empty() || elements.size() == 1) return;

    quick_sort_partition(elements, 0, elements.size() - 1);
}


// ------

/*
 * Merge sort -- O(N log N)
 * @param elements: a reference to vector of elements to be sorted
 *
 * Done in-place, so we return nothing
 */
template<typename T>
void merge(std::vector<T>& left, std::vector<T>& right, std::vector<T>& merged) {
    size_t left_size = left.size();
    size_t right_size = right.size();

    size_t i = 0;
    size_t j = 0;
    size_t k = 0;

    while (i < left_size && j < right_size) {
        if (left[i] < right[j]) {
            merged[k++] = left[i++];
        } else {
            merged[k++] = right[j++];
        }
    }

    while (i < left_size) {
        merged[k++] = left[i++];
    }

    while (j < right_size) {
        merged[k++] = right[j++];
    }
}

template<typename T>
void merge_sort(std::vector<T>& elements) { // NOLINT
    if (elements.empty() || elements.size() == 1) return;

    size_t size = elements.size();
    size_t mid = size / 2;

    std::vector<T> left(mid);
    std::vector<T> right(size - mid);

    for (size_t i = 0; i < mid; ++i) {
        left[i] = elements[i];
    }

    for (size_t i = mid; i < size; ++i) {
        right[i - mid] = elements[i];
    }

    merge_sort(left);
    merge_sort(right);

    merge(left, right, elements);
}

// --------

/*
 * Test cases generator
 */
std::vector<int> unordered_vec_of_ints() {
    return std::vector<int>{4, 2, 5, 1, 3, 3, 6};
}

int main() {
    // due to the in-order approach we kinda repeat ourselves in the tests
    std::vector<int> numbers;

    std::cout << "All should print 1 i.e True" << std::endl;

    /* insertion sort */
    numbers = unordered_vec_of_ints();
    insertion_sort(numbers);
    std::cout << std::is_sorted(numbers.begin(), numbers.end()) << std::endl;

    /* bubble sort */
    numbers = unordered_vec_of_ints();
    bubble_sort(numbers);
    std::cout << std::is_sorted(numbers.begin(), numbers.end()) << std::endl;


    /* selection sort */
    numbers = unordered_vec_of_ints();
    selection_sort(numbers);
    std::cout << std::is_sorted(numbers.begin(), numbers.end()) << std::endl;

    /* naive quick sort */
    numbers = unordered_vec_of_ints();
    auto sorted = naive_quick_sort(numbers);
    std::cout << std::is_sorted(sorted.begin(), sorted.end()) << std::endl;

    /* quick sort */
    numbers = unordered_vec_of_ints();
    quick_sort(numbers);
    std::cout << std::is_sorted(numbers.begin(), numbers.end()) << std::endl;

    /* merge sort */
    numbers = unordered_vec_of_ints();
    merge_sort(numbers);
    std::cout << std::is_sorted(numbers.begin(), numbers.end()) << std::endl;

    std::cout << "All done!" << std::endl;

    return 0;
}
