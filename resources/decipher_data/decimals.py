import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

# Define shifts to try
def shifts(v):
    return [
        v,
        v * 10,
        v / 10 if v != 0 else 0,
        v * 100,
        v / 100 if v != 0 else 0
    ]

# Function to optimize one row (greedy decimal fix)
def fix_row_decimal(row_values):
    values = row_values.copy()
    best_values = values.copy()
    current_sum = sum(values)

    for i in range(len(values)):
        best_local = values[i]
        best_local_sum = abs(current_sum - 100)

        for option in shifts(values[i]):
            temp_values = best_values.copy()
            temp_values[i] = option
            temp_sum = temp_values.sum()
            diff = abs(temp_sum - 100)

            if diff < best_local_sum:
                best_local = option
                best_local_sum = diff

        best_values[i] = best_local

    return list(best_values) + [sum(best_values), abs(sum(best_values) - 100)]

# Parallel processing wrapper
def parallel_fix(data):
    with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        futures = [executor.submit(fix_row_decimal, row.values) for _, row in data.iterrows()]
        results = [f.result() for f in as_completed(futures)]
    return results

# Main execution block — THIS is the fix
if __name__ == '__main__':
    # Load data
    df = pd.read_csv("data.csv")

    # Numeric columns only
    numeric_cols = df.select_dtypes(include='number').columns.tolist()
    df_numeric = df[numeric_cols]

    print("Starting parallel decimal fix...")

    results = parallel_fix(df_numeric)

    # Build result DataFrame
    corrected_array = np.array(results)
    corrected_df = pd.DataFrame(corrected_array[:, :-2], columns=numeric_cols)
    corrected_df["row_sum"] = corrected_array[:, -2]
    corrected_df["diff_from_100"] = corrected_array[:, -1]

    # Save results
    corrected_df.to_csv("greedy_decimal_fix_parallel.csv", index=False)
    print("✅ Done! File saved as 'greedy_decimal_fix_parallel.csv'")
