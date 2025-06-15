import tkinter as tk #biblioteca para interfaces gráficas, ttk fornece widgets 
from tkinter import ttk, messagebox #messagebox, para mostrar mensagens de erro/sucesso
import matplotlib.pyplot as plt #para criar gráficos, como o de barras
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg #integra os gráficos do matplotlib no

from sequence import Sequence #armazenar e manipular as sequências
from lcsfinder import LCSFinder #calcular o LCS
from sequencealignment import SequenceAlignment #alinhar as sequências
from sequencedatabase import SequenceDataBase #armazenar e manipular as sequências

class GeneLCSApp(tk.Tk): #define a janela principal da aplicação
    def __init__(self,filename):
        super().__init__()
        self.title("Comparador de Genes LCS")
        self.geometry("1000x800")
        self.minsize(800,600)
        self.config(bg="#FFC0CB")

        self.database = SequenceDataBase() #onde são guardadas as sequências
        self.sequences_file = filename

        self.font_title = ("Comic Sans MS", 20, "bold")
        self.font_label = ("Comic Sans MS", 14)
        self.font_button = ("Comic Sans MS", 12, "bold")
        
        self.selected_sequences = [None, None, None] #até 3 genes a comparar
        self.sequence_alignment_obj = None #resultado do alinhamento
        self.lcs_seq = "" #a sequência LCS calculada

        self.init_main_screen() #mostra o ecrã inicial
        self.load_initial_sequences() #carrega o ficheiro FASTA inicial

    def load_initial_sequences(self): #tenta carregar o ficheiro FASTA passado ao iniciar o programa
        try:
            self.database.load_from_fasta(self.sequences_file)
            self.update_gene_dropdowns()
            messagebox.showinfo("Sucesso", f"Sequências carregadas de '{self.sequences_file}' com sucesso!")
        except Exception as e:
            messagebox.showerror("Erro", f"Erro ao carregar as sequências: {e}")

    def update_gene_dropdowns(self): #atualiza as opções nos menus suspensos (Comboboxes) com os IDS e descrições das sequências
        gene_options = [(seq.id, seq.description) for seq in self.database.sequences]
        formatted_options = [f"ID: {id} - {desc}" for id, desc in gene_options]
        formatted_options.insert(0, "")

        for combo in self.gene_combos:
            combo['values'] = formatted_options
            combo.set('')

        self.selected_sequences = [None, None, None]
        self.compare_button.config(state=tk.DISABLED)

    def init_main_screen(self):
        for widget in self.winfo_children(): #limpa o conteúdo da janela, destrói os widgets
            widget.destroy()

        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

        outer_frame = tk.Frame(self, bg="#FFC0CB")
        outer_frame.grid(row=0, column=0, sticky="nsew")
        outer_frame.grid_rowconfigure(0, weight=1)
        outer_frame.grid_columnconfigure(0, weight=1)

        frame = tk.Frame(outer_frame, bg="#FFC0CB", padx=20, pady=20)
        frame.grid(row=0, column=0)
        
        title = tk.Label(frame, text="Selecionar Genes para Comparação", font=self.font_title, bg="#FFC0CB")
        title.grid(row=0, column=0, columnspan=2, pady=10)
        title.configure(anchor="center")

        self.gene_combos = []

        for i in range(3):
            label = tk.Label(frame, text=f"Gene {i+1}:", font=self.font_label, bg="#FFC0CB", anchor="center")
            label.grid(row=i+1, column=0, sticky="e", padx=(10, 5), pady=5)

            combo = ttk.Combobox(frame, state="readonly", width=60)
            combo.grid(row=i+1, column=1, sticky="w", padx=(0, 10))
            combo.bind("<<ComboboxSelected>>", self.on_gene_selected)
            self.gene_combos.append(combo)

        frame.grid_columnconfigure(0, weight=1)
        frame.grid_columnconfigure(1, weight=1)

        import_button = tk.Button(frame, text="Importar FASTA Manualmente", font=self.font_button,
                                bg="#FF69B4", fg="white",
                                command=self.import_fasta_file)
        import_button.grid(row=4, column=0, columnspan=2, pady=10)

        self.compare_button = tk.Button(frame, text="Compute LCS", font=self.font_button,
                                        bg="#FF69B4", fg="white", state=tk.DISABLED,
                                        command=self.compute_lcs)
        self.compare_button.grid(row=5, column=0, columnspan=2, pady=20)

        exit_button = tk.Button(frame, text="Sair da Aplicação", font=self.font_button,
                                bg="#FF69B4", fg="white", command=self.quit)
        exit_button.grid(row=6, column=0, columnspan=2, pady=10)

        self.update_gene_dropdowns()


    def import_fasta_file(self): #abre um seletor de ficheiros, carrega novas sequencias e se houver duplicados, adiciona sufixos
        from tkinter import filedialog

        filepath = filedialog.askopenfilename(
            title="Selecionar ficheiro FASTA",
            filetypes=[("FASTA files", ".fasta *.fa *.fna"), ("All files", ".*")]
        )
        if not filepath:
            return  # Cancelado

        try:
            new_db = SequenceDataBase()
            new_db.load_from_fasta(filepath)

            existing_ids = set(seq.id for seq in self.database.sequences)

            for seq in new_db.sequences:
                original_id = seq.id
                suffix = 1
                # Se já existe, renomear
                while seq.id in existing_ids:
                    seq.id = f"{original_id}_copy{suffix}"
                    suffix += 1
                self.database.sequences.append(seq)
                existing_ids.add(seq.id)

            self.update_gene_dropdowns()
            messagebox.showinfo("Sucesso", "Sequências importadas com sucesso!")

        except Exception as e:
            messagebox.showerror("Erro", f"Erro ao importar: {e}")



    def on_gene_selected(self, event): #quando um gene é escolhido no combobox, é obtido o ID
        #guarda as sequências selecionadas, ativa o botão "Compute LCS" apenas se forem selecionadas
        #2 ou 3 sequências diferentes
        for i, combo in enumerate(self.gene_combos):
            val = combo.get()
            if val:
                seq_id_str = val.split(" - ")[0].replace("ID: ", "")
                self.selected_sequences[i] = self.database.get_sequence_by_id(seq_id_str)
            else:
                self.selected_sequences[i] = None

        selected = [s for s in self.selected_sequences if s is not None]
        unique_ids = set(s.id for s in selected)
        if 2 <= len(selected) <= 3 and len(unique_ids) == len(selected):
            self.compare_button.config(state=tk.NORMAL)
        else:
            self.compare_button.config(state=tk.DISABLED)

    def compute_lcs(self): #usa a classe LCS para calcular a LCS
        active = [s for s in self.selected_sequences if s is not None]
        if len(active) < 2:
            messagebox.showerror("Erro", "Selecione pelo menos dois genes diferentes.")
            return

        try:
            if len(active) == 2:
                finder = LCSFinder(active[0], active[1])
            else:
                finder = LCSFinder(active[0], active[1], active[2])

            alignment = finder.compute_lcs()
            #guarda o resultado em self.sequence_aligment_obj e a LCS em self.lcs_seq
            self.sequence_alignment_obj = alignment
            self.lcs_seq = finder.lcs_seq

            self.show_results_screen()

        except Exception as e:
            messagebox.showerror("Erro", f"Erro ao calcular LCS: {e}")

    def show_results_screen(self): #mostra a sequência LCS, a identidade e 3 botoes
        for widget in self.winfo_children():
            widget.destroy()

        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

        frame = tk.Frame(self, bg="#FFC0CB", padx=20, pady=20)
        frame.grid(row=0, column=0, sticky="nsew")

        title = tk.Label(frame, text="Resultados do LCS", font=self.font_title, bg="#FFC0CB")
        title.pack(pady=10)

        lcs_label = tk.Label(frame, text="Sequência LCS:", font=self.font_label, bg="#FFC0CB")
        lcs_label.pack(anchor="center")

        lcs_text = tk.Text(frame, height=2, wrap=tk.WORD, font=("Courier New", 14), bg="white")
        lcs_text.pack(fill=tk.X, pady=5)
        lcs_text.insert(tk.END, self.lcs_seq)
        lcs_text.config(state=tk.DISABLED)

        identity_label = tk.Label(frame, text="Identity (%):", font=self.font_label, bg="#FFC0CB")
        identity_label.pack(anchor="center")

        identity_text = tk.Text(frame, height=1, font=("Helvetica", 14), bg="white")
        identity_text.pack(fill=tk.X, pady=5)
        identity_val = self.sequence_alignment_obj.identity()
        identity_text.insert(tk.END, f"{identity_val}%")
        identity_text.config(state=tk.DISABLED)

        btn_frame = tk.Frame(frame, bg="#FFC0CB")
        btn_frame.pack(pady=15)

        btn_color_map = tk.Button(btn_frame, text="Ver Gráfico de Cores", font=self.font_button, bg="#FF69B4", fg="white",
                                  command=self.show_color_map)
        btn_color_map.pack(side=tk.LEFT, padx=10)

        btn_lcs_length = tk.Button(btn_frame, text="Ver Gráfico do Comprimento da LCS", font=self.font_button, bg="#FF69B4", fg="white",
                                   command=self.show_lcs_length_chart)
        btn_lcs_length.pack(side=tk.LEFT, padx=10)

        btn_aligned_seq = tk.Button(btn_frame, text="Ver Sequências Alinhadas", font=self.font_button, bg="#FF69B4", fg="white",
                                    command=self.show_aligned_sequences)
        btn_aligned_seq.pack(side=tk.LEFT, padx=10)

        self.display_frame = tk.Frame(frame, bg="#FFC0CB")
        self.display_frame.pack(fill=tk.BOTH, expand=True)

        back_btn = tk.Button(frame, text="Voltar ao Menu Inicial", font=self.font_button, bg="#FF69B4", fg="white",
                             command=self.init_main_screen)
        back_btn.pack(pady=10)

    def show_color_map(self): #mostra cada caracter das sequências alinhadas com cores
        for widget in self.display_frame.winfo_children():
            widget.destroy()

        # cria um frame para scroll horizontal
        container = tk.Frame(self.display_frame, bg="#FFC0CB")
        container.pack(fill=tk.BOTH, expand=True)

        # LEGENDA
        legend_frame = tk.Frame(container, bg="#FFC0CB")
        legend_frame.pack(side=tk.TOP, pady=5)

        # Caixa verde
        tk.Label(legend_frame, bg="green", width=2, height=1).pack(side=tk.LEFT, padx=2)
        tk.Label(legend_frame, text="Match", bg="#FFC0CB").pack(side=tk.LEFT, padx=5)
        # Caixa amarela
        tk.Label(legend_frame, bg="yellow", width=2, height=1).pack(side=tk.LEFT, padx=2)
        tk.Label(legend_frame, text="Missmatch", bg="#FFC0CB").pack(side=tk.LEFT, padx=5)
        # Caixa vermelha
        tk.Label(legend_frame, bg="orange", width=2, height=1).pack(side=tk.LEFT, padx=2)
        tk.Label(legend_frame, text="Gaps \"-\"", bg="#FFC0CB").pack(side=tk.LEFT, padx=5)

        canvas = tk.Canvas(container, bg="#FFC0CB", height=100)
        canvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        scrollbar_x = tk.Scrollbar(container, orient="horizontal", command=canvas.xview)
        scrollbar_x.pack(side=tk.BOTTOM, fill=tk.X)

        canvas.configure(xscrollcommand=scrollbar_x.set)

        color_map_frame = tk.Frame(canvas, bg="#FFC0CB")
        canvas.create_window((0, 0), window=color_map_frame, anchor="nw")

        seqs = [self.sequence_alignment_obj.aligned_seq1, self.sequence_alignment_obj.aligned_seq2]
        if self.sequence_alignment_obj.seq3:
            seqs.append(self.sequence_alignment_obj.aligned_seq3)

        alignment_len = len(seqs[0])
        # Nova lógica para colorir cada elemento individualmente
        for row, seq in enumerate(seqs):
            for col in range(alignment_len):
                char = seq[col]
                column_chars = [s[col] for s in seqs]
                color = "white"  # default

                if char == "-":
                    color = "orange" # cor para gaps
                else:
                    # Conta quantas vezes esta letra aparece na coluna (ignorando '-')
                    if char != "-":
                        count = column_chars.count(char)
                        if count > 1:
                            # Só as letras repetidas ficam verdes
                            color = "green"
                        else:
                            color = "yellow"

                lbl = tk.Label(color_map_frame, text=char, bg=color, fg="black", font=("Courier New", 10), width=2)
                lbl.grid(row=row, column=col, padx=0, pady=0)

        def on_configure(event):
            canvas.configure(scrollregion=canvas.bbox("all"))

        color_map_frame.bind("<Configure>", on_configure)

    def show_lcs_length_chart(self): #mostra gráfico de barras, comprimento LCS (barra verde)
        #comprimento da maior sequência original (barra rosa)
        for widget in self.display_frame.winfo_children():
            widget.destroy()

        # Recolher todas as sequências utilizadas no alinhamento
        original_sequences = []
        if self.sequence_alignment_obj.seq1:
            original_sequences.append(self.sequence_alignment_obj.seq1)
        if self.sequence_alignment_obj.seq2:
            original_sequences.append(self.sequence_alignment_obj.seq2)
        if self.sequence_alignment_obj.seq3:
            original_sequences.append(self.sequence_alignment_obj.seq3)

        lcs_length = len(self.lcs_seq)
        max_seq_length = max(len(seq.seq) for seq in original_sequences)


        fig, ax = plt.subplots(figsize=(6, 4))

        labels = ['Comparação']
        lcs_lengths = [lcs_length]
        max_lengths = [max_seq_length]

        x = range(len(labels))
        width = 0.35

        bars1 = ax.bar([i - width/2 for i in x], lcs_lengths, width, label="LCS", color="green")
        bars2 = ax.bar([i + width/2 for i in x], max_lengths, width, label="Maior sequência", color="pink")

        ax.set_ylabel('Comprimento')
        ax.set_title('LCS vs Maior Sequência')
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.legend()

        for bar in bars1 + bars2:
            height = bar.get_height()
            ax.annotate(f'{height}',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=8)

        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=self.display_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)



    def show_aligned_sequences(self): #mostra as sequências alinhadas dentro de caixas
        for widget in self.display_frame.winfo_children():
            widget.destroy()

        #frame container com altura fixa e scroll vertical
        container_height = 400 
        container = tk.Frame(self.display_frame, bg="#FFC0CB", height=container_height)
        container.pack(fill="both", expand=False)

        canvas = tk.Canvas(container, bg="#FFC0CB", height=container_height)
        scrollbar_y = tk.Scrollbar(container, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar_y.set)

        scrollbar_y.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)

        scrollable_frame = tk.Frame(canvas, bg="#FFC0CB")
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")

        def on_frame_configure(event):
            canvas.configure(scrollregion=canvas.bbox("all"))

        scrollable_frame.bind("<Configure>", on_frame_configure)

        aligned_seqs = [self.sequence_alignment_obj.aligned_seq1,
                        self.sequence_alignment_obj.aligned_seq2]
        if self.sequence_alignment_obj.seq3:
            aligned_seqs.append(self.sequence_alignment_obj.aligned_seq3)

        for i, seq in enumerate(aligned_seqs):
            box = tk.LabelFrame(scrollable_frame, text=f"Sequência alinhada {i+1}",
                                font=self.font_label, bg="#FFC0CB", fg="black", labelanchor="n")
            box.pack(fill=tk.X, padx=10, pady=10)

            text = tk.Text(box, height=5, font=("Courier New", 14), wrap=tk.NONE, bg="white")
            text.pack(fill="both", expand=True, padx=5, pady=5)
            text.insert(tk.END, seq)
            text.config(state=tk.DISABLED)

        back_btn = tk.Button(self.display_frame, text="Voltar ao Menu Inicial",
                     font=self.font_button, bg="#FF69B4", fg="white",
                     command=self.init_main_screen)
        back_btn.pack(pady=10, padx=10) 


if __name__ == "__main__": #inicializa a aplicação carregando o ficheiro 
    app = GeneLCSApp("sequences.fasta")
    app.mainloop()