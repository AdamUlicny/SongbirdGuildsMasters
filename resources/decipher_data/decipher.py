import pandas as pd
import itertools

# Load your CSV file
df = pd.read_csv("data.csv")

# Select numeric columns to test (can also specify manually if needed)
columns_to_test = df.select_dtypes(include='number').columns.tolist()

# Set maximum number of columns to remove in any combination
MAX_REMOVE = 10

best_subset = None
best_avg_diff = float('inf')
best_remaining_sum = None

# Generate all combinations of columns to remove (from 0 up to MAX_REMOVE)
for r in range(0, MAX_REMOVE + 1):
    for remove_combo in itertools.combinations(columns_to_test, r):
        # Remaining columns after removing this combo
        remaining_cols = [col for col in columns_to_test if col not in remove_combo]

        # Calculate row-wise sum of remaining columns
        row_sums = df[remaining_cols].sum(axis=1)

        # Calculate average absolute difference from 100
        avg_diff = (row_sums - 100).abs().mean()

        # Track best-performing combo
        if avg_diff < best_avg_diff:
            best_avg_diff = avg_diff
            best_subset = remove_combo
            best_remaining_sum = row_sums

# Output the best result
print("Best columns to remove:", best_subset)
print("Average deviation from 100 after removal:", round(best_avg_diff, 3))

# Save result with the updated row sums to CSV
df_result = df.copy()
df_result["remaining_sum"] = best_remaining_sum
df_result.to_csv("data_with_best_removal.csv", index=False)

# Save summary
summary = pd.DataFrame({
    "columns_removed": [best_subset],
    "average_diff_from_100": [best_avg_diff]
})
summary.to_csv("best_column_removal_summary.csv", index=False)

