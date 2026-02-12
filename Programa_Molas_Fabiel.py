import math
import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib.pyplot as plt

# =============================================================================
# BANCO DE DADOS (TABELAS)
# =============================================================================

# Tabela 4 - Parâmetros dos Materiais
# Estrutura: 'Nome': [(d_min, d_max, G(GPa), m, A(MPa), Custo, Multiplicador_Ssy)]
MATERIAIS_DB = {
    "Fio musical A228": [
        (0.0, 0.8, 82.7, 0.145, 2211, 2.6, 0.45),
        (0.8, 1.61, 81.7, 0.145, 2211, 2.6, 0.45),
        (1.61, 3.0, 81.0, 0.145, 2211, 2.6, 0.45),
        (3.0, 6.5, 80.0, 0.145, 2211, 2.6, 0.45)
    ],
    "Fio duro estirado A227": [
        (0.7, 0.8, 80.7, 0.19, 1783, 1.0, 0.45),
        (0.8, 1.6, 80.0, 0.19, 1783, 1.0, 0.45),
        (1.6, 3.0, 79.3, 0.19, 1783, 1.0, 0.45),
        (3.0, 12.7, 78.6, 0.19, 1783, 1.0, 0.45)
    ],
    "Revenido em óleo A229": [
        (0.5, 12.7, 77.2, 0.187, 1855, 1.3, 0.50)
    ],
    "Cromo-vanádio A232": [
        (0.8, 11.1, 77.2, 0.168, 2005, 3.1, 0.50)
    ],
    "Cromo-silício A401": [
        (1.6, 9.5, 77.2, 0.108, 1974, 4.0, 0.50)
    ],
    "Inoxidável A313": [
        (0.3, 2.5, 69.0, 0.146, 1867, 7.6, 0.35), # Custo médio 7.6-11
        (2.5, 5.0, 69.0, 0.263, 2065, 7.6, 0.35),
        (5.0, 10.0, 69.0, 0.478, 2911, 7.6, 0.35)
    ],
    "Fósforo-bronze B159": [
        (0.1, 0.6, 41.4, 0.0, 1000, 8.0, 0.35),
        (0.6, 2.0, 41.4, 0.028, 913, 8.0, 0.35),
        (2.0, 7.5, 41.4, 0.064, 932, 8.0, 0.35)
    ]
}

# Tabela 1 - Tipos de Extremidade
TIPOS_EXTREMIDADE = [
    "Simples",
    "Simples e retificada",
    "Em esquadro ou fechada",
    "Em esquadro e retificada"
]

# Tabela 2 - Condição de Extremidade (Valor de alpha)
CONDICOES_EXTREMIDADE = {
    "Superfícies planas paralelas (fixas)": 0.5,
    "Uma fixa, uma pivotada": 0.707,
    "Ambas pivotadas": 1.0,
    "Uma engastada, uma livre": 2.0
}

# Tabela 3 - Tipo de Tratamento de Fadiga (Ssa, Ssm)
JATEAMENTO = {
    "Sem jateamento": {"Ssa": 241, "Ssm": 534},
    "Com jateamento": {"Ssa": 398, "Ssm": 379}
}

# Peso específico aproximado do aço
GAMMA_ACO = 7.7e-5  # N/mm^3 (aprox para aço)

# =============================================================================
# CÁLCULOS (Classe MolaCalculada)
# =============================================================================
class MolaCalculada:
    def __init__(self, d, D, material_nome, dados_entrada):
        self.d = d
        self.D = D
        self.C = D / d
        self.material_nome = material_nome
        self.inputs = dados_entrada
        self.valida = True
        self.mensagens_erro = []
        self.calculo_fadiga_ativo = False

        # Propriedades do material baseadas no diâmetro
        self.props = self._get_propriedades_material()

        if not self.props:
            self.valida = False
            self.mensagens_erro.append(f"Diâmetro {d}mm fora da faixa para {material_nome}")
            return

        self.G = self.props['G'] * 1000
        self.A = self.props['A']
        self.m = self.props['m']
        self.custo = self.props['custo']
        self.mult_Ssy = self.props['mult_ssy']

        # Cálculos Geométricos e de Força
        self._calcular_geometria()
        self._calcular_tensoes_e_seguranca()
        self._calcular_frequencia()
        self._calcular_fom()
        self._validar()

    # Material_nome é buscado da Tabela 1
    def _get_propriedades_material(self):
        ranges = MATERIAIS_DB.get(self.material_nome, [])
        for (dmin, dmax, G, m, A, custo, mult) in ranges:
            # Verifica faixa de "d" da tabela
            if dmin <= self.d <= dmax + 0.001:
                return {'G': G, 'm': m, 'A': A, 'custo': custo, 'mult_ssy': mult}
        return None

    def _calcular_geometria(self):
        # Fator Bergsträsser
        self.KB = (4 * self.C + 2) / (4 * self.C - 3)

        # Resistências
        self.Sut = self.A / (self.d ** self.m)
        self.Ssy = self.mult_Ssy * self.Sut

        # Definição de Na e Nt baseada na Tabela 1
        NT = self.inputs['NT']
        self.NT = NT
        tipo_ext = self.inputs['tipo_ext']

        if tipo_ext == "Simples":
            self.Na = NT
            self.Ne = 0
        elif tipo_ext == "Simples e retificada":
            self.Na = NT - 1
            self.Ne = 1
        elif tipo_ext == "Em esquadro ou fechada":
            self.Na = NT - 2
            self.Ne = 2
        else: # Em esquadro e retificada
            self.Na = NT - 2
            self.Ne = 2

        # Razão de mola k
        self.k = (self.d**4 * self.G) / (8 * self.D**3 * self.Na)

        # Cálculo de Ls e Passo
        self.L0 = self.inputs['L0']

        if tipo_ext == "Simples":
            self.rho = (self.L0 - self.d) / self.Na
            self.Ls = self.d * (NT + 1)
        elif tipo_ext == "Simples e retificada":
            self.rho = self.L0 / (self.Na + 1)
            self.Ls = self.d * NT
        elif tipo_ext == "Em esquadro ou fechada":
            self.rho = (self.L0 - 3 * self.d) / self.Na
            self.Ls = self.d * (NT + 1)
        else: # Em esquadro e retificada
            self.rho = (self.L0 - 2 * self.d) / self.Na
            self.Ls = self.d * NT

        # Deflexão até o sólido
        self.ys = self.L0 - self.Ls

    def _calcular_tensoes_e_seguranca(self):
        Fmax = self.inputs['Fmax']
        Fmin = self.inputs['Fmin']

        # --- LÓGICA DE FADIGA ---
        if Fmax is not None and Fmin is not None:
            self.calculo_fadiga_ativo = True
            # Forças e Tensões Alternada e Média
            self.Fa = (Fmax - Fmin) / 2
            self.Fm = (Fmax + Fmin) / 2
            self.tau_a = (self.KB * 8 * self.Fa * self.D) / (math.pi * self.d**3)
            self.tau_m = (self.KB * 8 * self.Fm * self.D) / (math.pi * self.d**3)

            # Tensão e Força no Fechamento Sólido
            xi = 0.15
            self.Fs = (1 + xi) * Fmax
            self.tau_s = (self.KB * 8 * self.Fs * self.D) / (math.pi * self.d**3)

            # Da tabela 3
            f_dados = JATEAMENTO[self.inputs['jateamento']]
            self.Ssa = f_dados['Ssa']
            self.nf = self.Ssa / self.tau_a if self.tau_a > 0 else 999.0
            self.ns = self.Ssy / self.tau_s
        else:
            # --- CÁLCULO ESTÁTICO (SEM FADIGA) ---
            self.calculo_fadiga_ativo = False
            self.Fs = self.k * self.ys # Força no fechamento sólido
            self.tau_s = (self.KB * 8 * self.Fs * self.D) / (math.pi * self.d**3)

            self.ns = self.Ssy / self.tau_s if self.tau_s > 0 else 0
            self.nf = 999.0 # Infinito
            self.tau_a = 0
            self.tau_m = 0

    def _calcular_frequencia(self):
        # Peso da mola (Gamma é omitido para aço, usando constante global já definida)
        self.W = GAMMA_ACO * (math.pi**2 * self.d**2 * self.NT * self.D) / 4
        # Frequência fundamental
        g = 9810
        self.f = 0.5 * math.sqrt((self.k * g) / self.W) if self.W > 0 else 0

    def _calcular_fom(self):
        # Figura de mérito
        self.fom = self.custo * self.W

    def _validar(self):
        # Validações

        # Índice de mola
        if not (4 <= self.C <= 12):
            self.valida = False
            self.mensagens_erro.append(f"C={self.C:.2f} fora de 4-12")

        # Espiras ativas
        if not (3 <= self.Na <= 15):
            self.valida = False
            self.mensagens_erro.append(f"Na={self.Na:.2f} fora de 3-15")

        # Flambagem
        alpha = CONDICOES_EXTREMIDADE[self.inputs['condicao_ext']]
        L0_crit = 2.63 * self.D / alpha
        if self.L0 >= L0_crit:
            self.valida = False
            self.mensagens_erro.append(f"L0 instável (Flambagem) > {L0_crit:.1f}")

        # Limitações de usuário
        if self.inputs.get('max_L0') and self.L0 > self.inputs['max_L0']:
            self.valida = False
            self.mensagens_erro.append("L0 excede limite do usuário")

        if self.inputs.get('max_Ls') and self.Ls > self.inputs['max_Ls']:
            self.valida = False
            self.mensagens_erro.append("Ls excede limite do usuário")

        # Fatores de Segurança
        if self.nf < 1.5:
            self.valida = False
            self.mensagens_erro.append(f"nf {self.nf:.2f} < 1.5")

        if self.ns < 1.2:
            self.valida = False
            self.mensagens_erro.append(f"ns {self.ns:.2f} < 1.2")


# =============================================================================
# INTERFACE GRÁFICA (AppMolas)
# =============================================================================
class AppMolas:
    def __init__(self, root):
        self.root = root
        self.root.title("Cálculo de Molas Helicoidais")

        # Frame Principal
        frame = ttk.Frame(root, padding="15")
        frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Inputs Gerais
        ttk.Label(frame, text="Dados Geométricos", font=('Arial', 10, 'bold')).grid(row=0, column=0, sticky=tk.W, pady=5)

        ttk.Label(frame, text="Diâmetros fio (d) [mm]:").grid(row=1, column=0, sticky=tk.W)
        self.entry_d = ttk.Entry(frame, width=30)
        self.entry_d.insert(0, "0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0")
        self.entry_d.grid(row=1, column=1)

        ttk.Label(frame, text="Diâmetro Médio (D) [mm]:").grid(row=2, column=0, sticky=tk.W)
        self.entry_D = ttk.Entry(frame); self.entry_D.grid(row=2, column=1)

        ttk.Label(frame, text="Total Espiras (NT):").grid(row=3, column=0, sticky=tk.W)
        self.entry_NT = ttk.Entry(frame); self.entry_NT.grid(row=3, column=1)

        ttk.Label(frame, text="Comp. Livre (L0) [mm]:").grid(row=4, column=0, sticky=tk.W)
        self.entry_L0 = ttk.Entry(frame); self.entry_L0.grid(row=4, column=1)

        # Fadiga (Checkbutton)
        ttk.Separator(frame, orient='horizontal').grid(row=5, column=0, columnspan=2, sticky='ew', pady=10)

        self.var_fadiga = tk.BooleanVar(value=False)
        self.chk_fadiga = ttk.Checkbutton(frame, text="Deseja realizar análise de Fadiga?",
                                          variable=self.var_fadiga, command=self.toggle_fadiga)
        self.chk_fadiga.grid(row=6, column=0, columnspan=2, sticky=tk.W)

        # Container para inputs de fadiga (Fmax, Fmin, Jateamento)
        self.widgets_fadiga = []

        lbl_fmax = ttk.Label(frame, text="Força Máxima (Fmax) [N]:")
        ent_fmax = ttk.Entry(frame)
        self.entry_Fmax = ent_fmax
        self.widgets_fadiga.extend([(lbl_fmax, 7, 0), (ent_fmax, 7, 1)])

        lbl_fmin = ttk.Label(frame, text="Força Mínima (Fmin) [N]:")
        ent_fmin = ttk.Entry(frame)
        self.entry_Fmin = ent_fmin
        self.widgets_fadiga.extend([(lbl_fmin, 8, 0), (ent_fmin, 8, 1)])

        lbl_jat = ttk.Label(frame, text="Jateamento:")
        cmb_jat = ttk.Combobox(frame, values=list(JATEAMENTO.keys()), state="readonly")
        cmb_jat.current(0)
        self.combo_jat = cmb_jat
        self.widgets_fadiga.extend([(lbl_jat, 9, 0), (cmb_jat, 9, 1)])

        # Inputs Finais
        ttk.Separator(frame, orient='horizontal').grid(row=10, column=0, columnspan=2, sticky='ew', pady=10)

        ttk.Label(frame, text="Tipo de Extremidade:").grid(row=11, column=0, sticky=tk.W)
        self.combo_ext = ttk.Combobox(frame, values=TIPOS_EXTREMIDADE, state="readonly")
        self.combo_ext.current(3)
        self.combo_ext.grid(row=11, column=1)

        ttk.Label(frame, text="Condição de Extremidade:").grid(row=12, column=0, sticky=tk.W)
        self.combo_cond = ttk.Combobox(frame, values=list(CONDICOES_EXTREMIDADE.keys()), state="readonly")
        self.combo_cond.current(0)
        self.combo_cond.grid(row=12, column=1)

        ttk.Label(frame, text="Material:").grid(row=13, column=0, sticky=tk.W)
        self.combo_mat = ttk.Combobox(frame, values=["Automático"] + list(MATERIAIS_DB.keys()), state="readonly")
        self.combo_mat.current(0)
        self.combo_mat.grid(row=13, column=1)

        ttk.Button(frame, text="CALCULAR", command=self.calcular).grid(row=15, column=0, columnspan=2, pady=15)

        self.btn_grafico = ttk.Button(frame, text="Ver Gráfico de Segurança", command=self.plotar_grafico, state='disabled')
        self.btn_grafico.grid(row=16, column=0, columnspan=2, pady=5)

        self.text_result = tk.Text(frame, height=15, width=60)
        self.text_result.grid(row=16, column=0, columnspan=2)

        self.toggle_fadiga()

    def toggle_fadiga(self):
        ativo = self.var_fadiga.get()
        for widget, r, c in self.widgets_fadiga:
            if ativo:
                widget.grid(row=r, column=c, sticky=tk.W if c==0 else tk.E, pady=2)
            else:
                widget.grid_remove()

# =============================================================================
# Cérebro do Programa (inputs, salva, filtra, cria)
# =============================================================================
    def calcular(self):
        try:
            # Coleta de dados
            diametros = [float(x.strip()) for x in self.entry_d.get().split(',')]
            inputs = {
                'D': float(self.entry_D.get()),
                'NT': float(self.entry_NT.get()),
                'L0': float(self.entry_L0.get()),
                'tipo_ext': self.combo_ext.get(),
                'condicao_ext': self.combo_cond.get(),
                'max_L0': None, 'max_Ls': None
            }

            # Lógica Condicional de Fadiga
            if self.var_fadiga.get():
                # Se marcado, TENTA pegar os valores. Se estiver vazio, vai dar erro de conversão (intencional)
                val_fmax = self.entry_Fmax.get()
                val_fmin = self.entry_Fmin.get()

                if not val_fmax or not val_fmin:
                    messagebox.showwarning("Atenção", "Para análise de fadiga, Fmax e Fmin são obrigatórios.")
                    return

                inputs['Fmax'] = float(val_fmax)
                inputs['Fmin'] = float(val_fmin)
                inputs['jateamento'] = self.combo_jat.get()
            else:
                # Se desmarcado, envia None
                inputs['Fmax'] = None
                inputs['Fmin'] = None
                inputs['jateamento'] = "Sem jateamento" # Irrelevante

            material_selecionado = self.combo_mat.get()

            molas_calculadas = []

            # Loop principal de cálculo (Teste de todos os materiais caso não possua)
            materiais_para_testar = list(MATERIAIS_DB.keys()) if material_selecionado == "Automático" else [material_selecionado]

            # Cria cada mola
            for mat in materiais_para_testar:
                for d in diametros:
                    mola = MolaCalculada(d, inputs['D'], mat, inputs)

            # Coloca cada mola na lista das calculadas
                    molas_calculadas.append(mola)

            # Captura dados para o gráfico
            self.dados_grafico = {'d': [], 'ns': [], 'nf': [], 'modo_fadiga': False}

            # Ordena por diâmetro para o gráfico não ficar riscado
            molas_ordenadas = sorted(molas_calculadas, key=lambda x: x.d)

            for m in molas_ordenadas:
                self.dados_grafico['d'].append(m.d)
                self.dados_grafico['ns'].append(m.ns)
                # Limita visualmente o nf infinito em 10 para não estragar a escala
                nf_vis = m.nf if m.nf < 10 else 10
                self.dados_grafico['nf'].append(nf_vis)
                self.dados_grafico['modo_fadiga'] = m.calculo_fadiga_ativo

            # Habilita o botão se houver dados
            if molas_calculadas:
                self.btn_grafico['state'] = 'normal'

            # Filtra válidas para seleção da melhor
            molas_validas = [m for m in molas_calculadas if m.valida]

            self.text_result.delete(1.0, tk.END)

            if not molas_validas:
                self.text_result.insert(tk.END, "Nenhuma mola válida encontrada.\nERROS:\n")
                for m in molas_calculadas[:5]:
                    self.text_result.insert(tk.END, f"d={m.d}: {', '.join(m.mensagens_erro)}\n")
            else:
                # Seleciona melhor mola pelo FOM (Menor custo*peso)
                melhor_mola = sorted(molas_validas, key=lambda x: x.fom)[0]

                res = f"--- RESULTADO (Modo: {'Fadiga' if melhor_mola.calculo_fadiga_ativo else 'Estático'}) ---\n"
                res = f"--- MELHOR MOLA SELECIONADA ---\n"
                res += f"Material: {melhor_mola.material_nome}\n"
                res += f"Diâmetro fio (d): {melhor_mola.d} mm\n"
                res += f"Índice (C): {melhor_mola.C:.2f}\n"
                res += f"Espiras ativas (Na): {melhor_mola.Na:.2f}\n"
                res += f"Comprimento Sólido (Ls): {melhor_mola.Ls:.2f} mm\n"
                res += f"Razão (k): {melhor_mola.k:.2f} N/mm\n"

                if melhor_mola.calculo_fadiga_ativo:
                    res += f"Fadiga (nf): {melhor_mola.nf:.2f}\n"
                    res += f"Estático (ns): {melhor_mola.ns:.3f}\n"
                else:
                    res += f"Fadiga: NÃO CALCULADA\n"
                    res += f"Estático (ns): {melhor_mola.ns:.3f} (no fechamento sólido)\n"
                    res += f"Força para fechar (Fs): {melhor_mola.Fs:.1f} N\n"

                res += f"FOM: {melhor_mola.fom:.4f}\n"

                # Aviso de Frequência
                res += f"\nAVISO: Se freq. trabalho > {melhor_mola.f/20:.1f} Hz, redesenhar.\n"

                self.text_result.insert(tk.END, res)


        except ValueError:
            messagebox.showerror("Erro", "Verifique se todos os campos numéricos estão preenchidos corretamente.")

    def plotar_grafico(self):
            if not hasattr(self, 'dados_grafico'):
                return

            d = self.dados_grafico['d']
            ns = self.dados_grafico['ns']
            nf = self.dados_grafico['nf']
            modo_fadiga = self.dados_grafico['modo_fadiga']

            plt.figure(figsize=(10, 6))

            # Linha do Fator Estático (ns) - Azul
            plt.plot(d, ns, marker='o', label='Fator Estático (ns)', color='blue', linestyle='-')
            plt.axhline(y=1.2, color='blue', linestyle='--', alpha=0.5, label='Min Estático (1.2)')

            # Linha do Fator de Fadiga (nf) - Verde (Só se existir)
            if modo_fadiga:
                plt.plot(d, nf, marker='s', label='Fator Fadiga (nf)', color='green', linestyle='-')
                plt.axhline(y=1.5, color='green', linestyle='--', alpha=0.5, label='Min Fadiga (1.5)')

            # Área Vermelha (Perigo abaixo de 1.2)
            plt.axhspan(0, 1.2, color='red', alpha=0.1)
            plt.text(min(d), 0.5, "ZONA DE FALHA", color='red', fontsize=12, fontweight='bold')

            plt.title(f"Sensibilidade: Diâmetro do Fio vs. Segurança")
            plt.xlabel("Diâmetro do Fio (d) [mm]")
            plt.ylabel("Fator de Segurança")

            # Ajusta limite Y para ficar bonito (pega o maior valor entre ns e nf, mas trava em 6 se for muito alto)
            max_y = max(max(ns), max(nf) if modo_fadiga else 0)
            plt.ylim(0, min(max_y + 1, 8))

            plt.grid(True, which='both', linestyle='--', linewidth=0.5)
            plt.legend()
            plt.tight_layout()
            plt.show()

if __name__ == "__main__":
    root = tk.Tk()
    app = AppMolas(root)
    root.mainloop()
