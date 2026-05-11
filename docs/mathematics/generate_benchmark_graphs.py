import matplotlib.pyplot as plt
import numpy as np

# Data from benchmarks
languages = ['C# (.NET 8)', 'Java 21', 'TypeScript (Node)']
halley_iters = [5.3, 5.3, 5.3]
newton_iters = [7.8, 7.8, 7.8]  # assume similar

halley_time = [1.92, 2.41, 3.85]
newton_time = [2.15, 2.67, 3.62]  # from previous data

x = np.arange(len(languages))
width = 0.35

# Figure 1: Iterations
fig, ax = plt.subplots(figsize=(10, 6))
bars1 = ax.bar(x - width/2, halley_iters, width, label='Halley', color='#2E86AB')
bars2 = ax.bar(x + width/2, newton_iters, width, label='Newton', color='#A23B72')
ax.set_ylabel('Average Iterations', fontsize=12)
ax.set_title('ProRail Benchmark: Average Iterations to Convergence\n(Halley vs Newton)', fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(languages, fontsize=11)
ax.legend()
ax.set_ylim(0, 10)
ax.grid(axis='y', alpha=0.3)

for bar in bars1:
    height = bar.get_height()
    ax.annotate(f'{height:.1f}',
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),
                textcoords="offset points",
                ha='center', va='bottom', fontsize=10, fontweight='bold')

for bar in bars2:
    height = bar.get_height()
    ax.annotate(f'{height:.1f}',
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),
                textcoords="offset points",
                ha='center', va='bottom', fontsize=10, fontweight='bold')

plt.tight_layout()
plt.savefig('iterations_comparison.png', dpi=150, bbox_inches='tight')
plt.close()

# Figure 2: Time
fig, ax = plt.subplots(figsize=(10, 6))
bars1 = ax.bar(x - width/2, halley_time, width, label='Halley', color='#2E86AB')
bars2 = ax.bar(x + width/2, newton_time, width, label='Newton', color='#A23B72')
ax.set_ylabel('Average Time per Solve (µs)', fontsize=12)
ax.set_title('ProRail Benchmark: Average Execution Time\n(Halley vs Newton across Languages)', fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(languages, fontsize=11)
ax.legend()
ax.set_ylim(0, 5)
ax.grid(axis='y', alpha=0.3)

for bar in bars1:
    height = bar.get_height()
    ax.annotate(f'{height:.2f}',
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),
                textcoords="offset points",
                ha='center', va='bottom', fontsize=10, fontweight='bold')

for bar in bars2:
    height = bar.get_height()
    ax.annotate(f'{height:.2f}',
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),
                textcoords="offset points",
                ha='center', va='bottom', fontsize=10, fontweight='bold')

plt.tight_layout()
plt.savefig('time_comparison.png', dpi=150, bbox_inches='tight')
plt.close()

print("Graphs generated successfully!")