import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os

# Importar as Classes
from sequence import Sequence
from lcsfinder import LCSFinder
from sequencealignment import SequenceAlignment
from sequencedatabase import SequenceDataBase

class GeneLCSApp:
    def __init__(self, master, filename):
        self.master = master
        master.title("Comparador de Genes LCS")
        master.geometry("800x700")

        self.database = SequenceDataBase()
        self.selected_sequences = [None, None, None]

        # --- Frame Principal ---
        self.main_frame = tk.Frame(master, padx=10, pady=10)
        self.main_frame.pack(fill=tk.BOTH, expand=True)

        # --- Seleção de Genes ---
        self.select_genes_frame = tk.LabelFrame(self.main_frame, text="Selecionar Genes para Comparação", padx=10, pady=10)
        self.select_genes_frame.pack(fill=tk.X, pady=10)

        # Dropdowns para seleção
        self.gene_labels = []
        self.gene_combos = []
        for i in range(3):
            label = tk.Label(self.select_genes_frame, text=f"Gene {i+1}:")
            label.grid(row=i, column=0, padx=5, pady=5, sticky="w")
            self.gene_labels.append(label)

            combo = ttk.Combobox(self.select_genes_frame, state="readonly")
            combo.grid(row=i, column=1, padx=5, pady=5, sticky="ew")
            combo.bind("<<ComboboxSelected>>", self.on_gene_selected)
            self.gene_combos.append(combo)

        self.select_genes_frame.grid_columnconfigure(1, weight=1)

        # Botão de Comparar
        self.compare_button = tk.Button(self.main_frame, text="Comparar Genes (LCS)", command=self.compare_genes, state=tk.DISABLED)
        self.compare_button.pack(pady=10)

        # Resultados LCS
        self.results_frame = tk.LabelFrame(self.main_frame, text="Resultados do LCS", padx=10, pady=10)
        self.results_frame.pack(fill=tk.BOTH, expand=True, pady=10)

        self.lcs_label = tk.Label(self.results_frame, text="Sequência LCS:")
        self.lcs_label.pack(anchor="w")
        self.lcs_text = tk.Text(self.results_frame, height=2, wrap=tk.WORD)
        self.lcs_text.pack(fill=tk.X, pady=5)

        self.aligned_label = tk.Label(self.results_frame, text="Sequências Alinhadas:")
        self.aligned_label.pack(anchor="w")
        self.aligned_text = tk.Text(self.results_frame, height=8, wrap=tk.WORD)
        self.aligned_text.pack(fill=tk.BOTH, expand=True, pady=5)

        # Área para Gráficos
        self.chart_frame = tk.LabelFrame(self.main_frame, text="Visualização Gráfica", padx=10, pady=10)
        self.chart_frame.pack(fill=tk.BOTH, expand=True, pady=10)
        self.fig, self.ax = plt.subplots(figsize=(6, 4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.chart_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill=tk.BOTH, expand=True)

        # --- Carregar Sequências Automaticamente ---
        self.sequences_file = filename
        self.load_initial_sequences()

    def load_initial_sequences(self):
        if os.path.exists(self.sequences_file):
            try:
                self.database.load_from_fasta(self.sequences_file)
                self.update_gene_dropdowns()
                messagebox.showinfo("Sucesso", f"Sequências carregadas de '{self.sequences_file}' com sucesso!")
            except Exception as e:
                messagebox.showerror("Erro", f"Erro ao carregar as sequências de '{self.sequences_file}': {e}")
                self.clear_dropdowns()
        else:
            messagebox.showerror("Erro", f"O arquivo de sequências '{self.sequences_file}' não foi encontrado. Por favor, crie-o no mesmo diretório.")

    def update_gene_dropdowns(self):
        gene_options = [(seq.id, seq.description) for seq in self.database.sequences]
        formatted_options = [f"ID: {id} - {desc}" for id, desc in gene_options]

        for combo in self.gene_combos:
            combo['values'] = formatted_options
            combo.set('') # Limpa a seleção anterior

        self.selected_sequences = [None, None, None] # Reset selected sequences
        self.check_comparison_readiness()

    def on_gene_selected(self, event):
        for i, combo in enumerate(self.gene_combos):
            selected_value = combo.get()
            if selected_value:
                seq_id_str = selected_value.split(" - ")[0].replace("ID: ", "")
                self.selected_sequences[i] = self.database.get_sequence_by_id(seq_id_str)
            else:
                self.selected_sequences[i] = None
        self.check_comparison_readiness()

    def check_comparison_readiness(self):
        num_selected = sum(1 for s in self.selected_sequences if s is not None)
        selected_ids = set(s.id for s in self.selected_sequences if s is not None)

        if 2 <= num_selected <= 3 and len(selected_ids) == num_selected:
            self.compare_button.config(state=tk.NORMAL)
        else:
            self.compare_button.config(state=tk.DISABLED)

    def compare_genes(self):
        seq1 = self.selected_sequences[0]
        seq2 = self.selected_sequences[1]
        seq3 = self.selected_sequences[2]

        active_sequences = [s for s in [seq1, seq2, seq3] if s is not None]

        if len(active_sequences) != len(set(s.id for s in active_sequences)):
            messagebox.showerror("Erro de Seleção", "Por favor, selecione genes diferentes para comparação.")
            self.compare_button.config(state=tk.DISABLED)
            return

        if len(active_sequences) < 2:
            messagebox.showerror("Erro", "Por favor, selecione pelo menos duas sequências para comparar.")
            return

        try:
            if len(active_sequences) == 2:
                finder = LCSFinder(active_sequences[0], active_sequences[1])
            else: # len(active_sequences) == 3
                finder = LCSFinder(active_sequences[0], active_sequences[1], active_sequences[2])

            sequence_alignment_obj = finder.compute_lcs()

            self.lcs_text.delete("1.0", tk.END)
            self.lcs_text.insert(tk.END, finder.lcs_seq) 

            self.aligned_text.delete("1.0", tk.END)
            if sequence_alignment_obj:
                self.aligned_text.insert(tk.END, str(sequence_alignment_obj))

            self.plot_lcs_length(finder.lcs_seq)

            messagebox.showinfo("Comparação Concluída", "Cálculo do LCS realizado com sucesso!")

        except Exception as e:
            messagebox.showerror("Erro de Comparação", f"Ocorreu um erro ao comparar os genes: {e}")

    def plot_lcs_length(self, lcs_seq):
        self.ax.clear()
        length = len(lcs_seq)
        if length > 0:
            self.ax.bar(['LCS Length'], [length], color='skyblue')
            self.ax.set_ylabel('Length')
            self.ax.set_title('Length of Longest Common Subsequence (LCS)')
            self.ax.text(0, length, str(length), ha='center', va='bottom')
        else:
            self.ax.text(0.5, 0.5, "No LCS found or LCS length is 0", horizontalalignment='center', verticalalignment='center', transform=self.ax.transAxes)

        self.canvas.draw()

    def clear_dropdowns(self):
        for combo in self.gene_combos:
            combo['values'] = []
            combo.set('')
        self.selected_sequences = [None, None, None]
        self.compare_button.config(state=tk.DISABLED)
        self.lcs_text.delete("1.0", tk.END)
        self.aligned_text.delete("1.0", tk.END)
        self.ax.clear()
        self.canvas.draw()

if __name__ == "__main__":
    root = tk.Tk()
    app = GeneLCSApp(root)
    root.mainloop()


