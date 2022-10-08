#include "unity/unity.h"

#include "bruits.h"

void setUp(void)
{
}

void tearDown(void)
{
}

void test_minimum(void)
{
    TEST_ASSERT_EQUAL_UINT8(2, br_minimum(2, 4));
    TEST_ASSERT_EQUAL_FLOAT(2., br_minimum(2., 4.));
    TEST_ASSERT_EQUAL_FLOAT(9.1, br_minimum(9.1, 13));
}

void test_maximum(void)
{
    TEST_ASSERT_EQUAL_UINT8(4, br_maximum(2, 4));
    TEST_ASSERT_EQUAL_FLOAT(4., br_maximum(2., 4.));
    TEST_ASSERT_EQUAL_FLOAT(13, br_maximum(9.1, 13));
}

void test_clamp(void)
{
    TEST_ASSERT_EQUAL_UINT8(0, br_clamp(-23, 0, 1));
    TEST_ASSERT_EQUAL_UINT8(1, br_clamp(23, 0, 1));

    TEST_ASSERT_EQUAL_FLOAT(0.1, br_clamp(-23, 0.1, 23.));
    TEST_ASSERT_EQUAL_FLOAT(3.14, br_clamp(3.14, 0.1, 23.));
    TEST_ASSERT_EQUAL_FLOAT(23, br_clamp(68, 0.1, 23.));
}

int main()
{
    UNITY_BEGIN();
    RUN_TEST(test_minimum);
    RUN_TEST(test_maximum);
    RUN_TEST(test_clamp);
    return UNITY_END();
}
