# AeroOpt v0.4.0 — verificação da integração Alula

## Conclusão e alcance

Os modos `inviscid` e `viscous_incompressible` foram integrados ao runner
existente por processos Node.js isolados. As traduções de geometria, sinais,
condições, coeficientes e estados foram verificadas contra chamadas diretas
reproduzíveis ao mesmo upstream fixado. Um caso de dois elementos viscoso
convergiu e conservou seus avisos. **Isso verifica a integração de software;
não estabelece validação experimental, precisão viscosa multielemento ou
aptidão para otimização aerodinâmica.**

Não foram alterados Hess–Smith, `PanelSolver`/sua identidade 0.3.1, runner,
supervisor, `SolverCase`, `OperatingConditions`, coeficientes ou schemas neutros.
O pacote passa a 0.4.0, com novo adaptador/bridge próprios. A v0.3.2 está na
base `db30d8eafbc181e67f8c877426069b30a250699b`; o ramo desta entrega é
`feature/alula-integration`. Não há merge para main.

| Evidência | O que estabelece | O que permanece aberto |
|---|---|---|
| Fontes/API/licença inspecionadas e fixadas | Componentes e bytes efetivamente executados identificados | Autenticação criptográfica do autor; licença global do AeroOpt anterior |
| Paridade com chamada direta upstream | Bridge conserva números, Cp, elementos e estados da API | Exatidão aerodinâmica do próprio upstream |
| Integração independente de pressões e transferência de origem | Sinais de CL/CM, origem e normalização aplicados corretamente | Sensibilidade geral a malha/forma e modelos físicos |
| Repetições e regressão | Consistência local e compatibilidade nos casos exercitados | Igualdade entre runtimes/SO ou todo o domínio de geometrias |
| Resíduos BL/esteira e admissibilidade | Aceitação numérica/modelo segundo gates declarados | Precisão experimental, interação de esteiras, separação forte e estol |

## Proveniência e instalação

Upstream [kennyjensen/alula](https://github.com/kennyjensen/alula), commit
**`160f43713a1c6eafb0a4b0f035f5b0bb70b016f8`**, de 2026-10-03. Foi inspecionado
o HEAD atual, suas APIs, retorno de cada modo, preparação/topologia dos
contornos, integração de forças, critérios de aceitação, dependências e avisos.
O snapshot isolado contém **47 arquivos / 510822 bytes**, incluindo **39
módulos JS** e oito arquivos de licença/metadados. Não contém UI, datasets,
PDFs ou módulos Euler. O lock com SHA-256 faz parte do pacote; nenhuma análise
efetua download. Instalação/limites/API/licenças estão em
[ALULA_INTEGRATION.md](ALULA_INTEGRATION.md).

Node exercitado: **v22.23.4 Windows x64**, portátil local em `.deps/`, obtido
do site oficial, ZIP e binário fixados por SHA-256. Python **3.13.15** na
`.venv` existente, Windows 11. Foram preservados os originais, checksums e
proveniência em `outputs/upstream_sources_v040/`. A instalação reproduzível
pode usar downloads explícitos ou uma cópia local conferida byte a byte.
O caminho `.deps/` não é versionado. Nenhuma dependência Python de execução
foi adicionada; não foram usados npm, navegador, WASM ou instalação global.

Alula/XFOIL são GPL-2.0-or-later, com avisos de Kenny Jensen e Mark Drela
preservados. Os wrappers JS desta integração têm esse identificador e COPYING
integral; não se atribuiu outra licença às fontes upstream nem ao restante
do AeroOpt. Uso não comercial não dispensa as condições da GPL. Referências
externas/datasets citados pelo upstream não foram incluídos. Não se afirma
que a separação em processos determine, por si só, o alcance jurídico de uma
distribuição combinada.

## Casos e convenções

Todos os perfis são NACA analíticos explicitamente discretizados; não há
coordenadas externas importadas ou ajustadas. Os dois métodos de referência
integram o mesmo problema/modelo upstream, sem usar o kernel Hess–Smith.

| Caso | Geometria e condições | Caminho exercitado |
|---|---|---|
| A | NACA 0012, 80 painéis, fuga afiada, alpha=0°, Re=1e6 como metadado, Mach=0 | JSON → bridge → `solveInviscid` |
| B | Mesmo perfil, alpha=4° | Bridge, três processos independentes e baselines diretos |
| C | NACA 0012 principal + NACA 4412 flap, 80 painéis/corpo, razão .3, LE flap (1.02,-.04), deflexão 25° horária, alpha=6°, Re=1e6 metadado, Mach=0 | `AlulaSolver` + `SolverRunner`, três tentativas + baselines fora das tentativas |
| D | NACA 0012, 80 painéis molhados, fuga finita padrão, alpha=0°, Re=1e6, Mach=0, Ncrit=9, trips=(1,1), 24 estações de esteira | Bridge → `solveViscousAssembly`, três repetições + baselines diretos |
| E | Perfil finito de D, alpha=2°, somente uma iteração Newton e zero atualizações de esteira | Não convergência deliberada, processo zero com coeficientes utilizáveis ausentes |
| V multielemento | Dois NACA 0012 afiados, 40 painéis/corpo, razão .3, LE flap (1.2,-.15), deflexão=0°, alpha=0°, Re=1e6, Mach=0, Ncrit=9, trips=(1,1) | Adaptador viscoso e runner reais |
| F finito | Montagem de C com endpoints de fuga distintos e tampa reta, streamfunction explícito | Adaptador invíscido, topologia de base finita |
| U | Montagem de C, Mach=.1 | Rejeição anterior ao processo numérico |

O runner continua sendo o contrato de **dois elementos**; A/B/D/E verificam o
bridge standalone, que suporta 1–6. A montagem Python não foi artificialmente
alterada para simular um único elemento. As coordenadas globais dos snapshots
são enviadas exatamente, principal antes do flap, CCW. Contornos horários são
rejeitados em vez de permitir a reordenação automática do upstream. Uma fuga
finita preserva os endpoints e acrescenta somente a tampa reta declarada.
Sua pressão de base é modelada pelo Alula; não foi qualificada como tratamento
experimental validado ou solução da limitação de fuga finita do Hess–Smith.

Todas as tabelas abaixo usam `c_main=1`, momento em `(0.25,0)`, **CCW positivo**,
normalizado por `q_inf*c_main²`. Alula retorna CM horário: compara-se o negativo
da chamada direta com o valor traduzido. Alpha define a velocidade livre
`(cos(alpha),sin(alpha))`. O teste com origem `(0.4,.03)` e sentido horário
confere a transferência de momento com as componentes de força integradas.
Não se trocam sinais empiricamente. As entradas físicas são explícitas.

## Resultados numéricos

Lote final: `outputs/alula_v040_final/`, execução supervisionada concluída em
**25.804 s**, sem incluir instalação. Números da integração, não medições
experimentais ou limites de erro físico:

| Caso | CL | CD utilizável | CM CCW | Resultado |
|---|---:|---:|---:|---|
| A | 3.50e-15 | ausente | -5.98e-16 | convergido/aceito |
| B | .476788529882641 | ausente | .006445853438383 | convergido/aceito |
| C | 3.119943444109189 | ausente | .806612864577950 | convergido/aceito, avisos de refinamento |
| D | 1.77e-15 | .005344683098715 | 2.61e-16 | convergido/aceito, aviso de refinamento |
| E | ausente | ausente | ausente | upstream `unconverged`; provisórios somente na saída bruta |
| V multielemento | .041225176883411 | .008513282903411 | .013006886843800 | convergido/aceito, avisos BL/refinamento |
| F finito | 3.162679516677851 | ausente | .813624263105625 | convergido/aceito, `ALULA_FINITE_BASE_RESEARCH` |
| U | ausente | ausente | ausente | `invalid_input`, Mach não reinterpretado |

As chamadas diretas reproduzíveis estão em
`examples/alula_upstream_baseline.mjs`. Elas não importam a tradução do bridge;
chamam as duas APIs upstream com os mesmos arrays/controles e deixam CM no
sentido horário original. Nos **dez pares** aceitos A + B×3 + C×3 + D×3, a
diferença de CL/CD/CM traduzidos foi zero em ponto flutuante. Os objetos completos
de elementos, incluindo coordenadas e distribuições Cp, também coincidiram
exatamente. Essa é uma verificação de integração: ambos executam as mesmas
equações upstream, logo a paridade não é uma validação aerodinâmica independente.

As constantes de regressão single-element vieram primeiro de execução direta
do upstream com seu gerador analítico `naca4('0012',80)`. As entradas AeroOpt
têm diferenças de arredondamento de ponto flutuante e deram diferença <=1e-14
nos coeficientes; os testes usam tolerância de integração 2e-11, documentada,
sem ajustar geometria. O caso C também coincide com a formulação LV independente
da v0.3.1 no mesmo nível de malha: CL=3.119943444109184 e CM=.806612864577952
no estudo anterior, diferença de poucos e-15. Não se inferiu disso precisão
física ou suficiência da malha: a v0.3.2 mostrou incrementos de refinamento
relevantes nessa montagem e candidatos próximos não resolvidos.

O caso D tem convergência simultânea de equações e esteira. No caso V, o máximo
resíduo das equações foi `6.217e-14`, da esteira `4.772e-7` (limites 1e-8/1e-6),
com sete estados de esteira e 62 iterações Newton acumuladas. O limite de 60
se aplica a cada solve interno, não à soma de todas as atualizações de esteira.
A maior espessura relativa foi `.0137381`, abaixo do limite amplo .1 do modelo.
Continuam os avisos upstream: malha com menos de 160 painéis por corpo,
separação local prevista e falta de referência viscosa multielemento independente.
Os valores são estimativas de pesquisa.

CD viscoso é obtido pelo déficit de momento da esteira. Em V, o diagnóstico
de arrasto por integral de pressão foi **-.00684470**, enquanto o CD de esteira
foi **+.00851328**. Não se substitui CD pelo diagnóstico de pressão, não se
zera arrasto negativo e não se soma duas estimativas diferentes. No invíscido,
`pressureDrag` permanece diagnóstico e CD continua ausente.

Em F, há 160 painéis molhados + duas tampas; os resíduos linear/Kutta ficaram
em `7.216e-16`/zero. A diferença CL pressão–circulação foi `.00126390` e o
arrasto integrado de pressão foi `-.00502603`. Isso não é evidência de arrasto
físico negativo nem validação da base. A limitação específica fica na saída,
assim como o aviso upstream de pressão de base modelada.

No caso E, o Node terminou com código zero, mas a solução ficou não convergida.
O bridge conservou coeficientes finitos da última iteração somente sob
`upstream.coefficients`, com rótulo provisório. O resíduo infinito de esteira
ficou `null` acompanhado de seu caminho/representação em
`nonfinite_diagnostics`. Os três coeficientes utilizáveis ficaram `null`.
Testes por injeção controlada também verificam falta de coeficiente/diagnóstico,
coeficiente não finito e estado `outside-model`: convergência matemática não
dispensa a rejeição por admissibilidade. Não se publicaram números inventados
como referência física nesses testes.

## Custos e consistência

Cinco processos supervisionados `node --version` levaram
`.2127, .1604, .1530, .1551, .1524 s`: mediana **.1551 s**. Essa é uma medida
de **startup + supervisor Windows**, não startup puro do motor V8. O JSON
também registra importação, integridade, tempo do solver, tempo de bridge e
uptime interno Node. A tabela usa medianas de três repetições:

| Caso | Solver upstream, s | Processo Node supervisionado, s | Runner Python total, s | Node menos solver, s |
|---|---:|---:|---:|---:|
| B, 80 painéis | .01891 | .23963 | — standalone | .22072 |
| D, 80 molhados + base | .66930 | .94108 | — standalone | .27178 |
| C, 160 painéis totais | .04460 | .42672 | 1.31853 | .38212 |

Em C, a diferença entre medianas runner/processo foi **.89181 s**, cobrindo
preparação, conferência de fontes/binário, validação, serialização, parsing e
inventário/auditoria Python. Subtrair também o baseline de startup dá resíduos
aproximados de overhead de `.06559 s` (B), `.11665 s` (D) e `.22699 s` (C).
Essas diferenças agregam importação, validação, bridge, JSON e I/O, e sofrem
ruído de medições; não são uma decomposição exata nem um benchmark puro de IPC.
O tempo do runner exclui a construção inicial do adaptador/sondagem de versão,
que pode ser reutilizado entre tentativas; cada solve ainda inicia outro Node.

V levou **4.477 s** no runner e F **1.275 s**. As amostras RSS do Node foram
aproximadamente **38–39 MiB** nos invíscidos, **54–55 MiB** em D e **70.6 MiB**
em V. São amostras antes da escrita final, não picos ou memória total da árvore;
excluem Python, supervisor, caches e eventual memória fora da amostragem.
Não houve otimização de worker persistente ou enfraquecimento de integridade.

Nas três repetições B/C/D, os spreads de CL/CD/CM presentes foram **zero** em
ponto flutuante. As três tentativas C tiveram o mesmo ID científico e UUIDs
distintos. Não se promete igualdade bit a bit em outro Node, V8, SO ou malha.
Configurações, limites e hashes estão conservados. Rodou-se a campanha sem
concorrência com a suíte de testes para não misturar seus custos.

## Regressão e artefatos

O modo padrão antes das alterações passou **403 testes**, ignorando 23 fontes
UIUC/Williams locais, em 86.44 s. Os **57 casos novos** passaram em 22.88 s.
A suíte padrão final passou **460 testes**, com os mesmos **23 skips**, sem
falhas, em **127.60 s**. A execução completa com fontes UIUC/Williams locais
passou **483 testes**, sem falhas ou skips, em **117.47 s**. Não houve warnings
pytest nas execuções finais; os avisos científicos dos solvers são preservados
nos resultados e não ocultados por esse resumo.

| Execução | Passaram | Skips | Falhas | Tempo, s |
|---|---:|---:|---:|---:|
| Base v0.3.2 | 403 | 23 | 0 | 86.44 |
| Somente novos testes Alula | 57 | 0 | 0 | 22.88 |
| Final padrão | 460 | 23 | 0 | 127.60 |
| Final com fontes locais | 483 | 0 | 0 | 117.47 |

Os novos testes cobrem:

- Round trip JSON, integração Python/Node e comparação upstream de um/dois
  elementos; resultados invíscidos e viscosos reais.
- Node ou snapshot ausente, fonte alterada, modos e parâmetros incompatíveis,
  limites de Re/alpha/painéis, Mach não zero, geometria e orientação inválidas.
- JSON inválido, coeficiente não finito/ausente, coeficiente provisório, falta
  de diagnósticos, modelo inadmissível e falsificação de flags/limitações.
- Identidade estável em diretórios diferentes, inputs congelados, hash de
  entrada e sem caminhos no ID; mudança de referência altera a identidade.
- Integração independente de pressão/forças/momento, conversão CW↔CCW e
  transferência de origem não padrão, sem alterar pontos.
- Caminhos com espaços/Unicode, Node absoluto, rejeição explícita de cwd longo,
  saída não zero e timeout com pai/descendente reais encerrados por Job Object.
- As regressões geométricas, de processos, Hess–Smith e verificação v0.3.2
  permanecem na suíte completa. Não se relaxaram suas tolerâncias.

Sem Node/snapshot opcional, os testes correspondentes dão skip explícito;
erros de contrato independentes dessas dependências continuam executáveis.
Uma instalação com fontes alteradas falha em vez de ser ignorada. Não se criou
teste que exigisse rede ou instalasse dependências durante pytest.

O lote final contém **228 artefatos** conferidos por SHA-256, seis snapshots de
código/lock/licença e seis tentativas do runner com inventários selados
conferidos. Não se acrescentaram baselines às pastas seladas: seus logs/resultados
ficam em diretórios separados. A igualdade dos dez conjuntos de elementos/Cp
foi conferida, além dos coeficientes. O inventário final tem SHA-256:

```text
ffffc3d27277f99f9b7497c5b03146cd3bcaa03271fb4b31f93794362b9ac37c
```

Também foi construído um wheel offline na `.venv` existente e conferida a
presença byte a byte de bridge/lock/COPYING. Ele não embute Node ou o snapshot
Alula. `pip check` passou, sem novas dependências Python de runtime.
Um resultado real da v0.3.2 também foi carregado pela API atual, preservando
status e identidade `aeroopt.hess_smith/0.3.1`.
O exemplo foi executado via `AttemptContext.run_process`, com todos os outputs dentro
do workspace. Logs de execução ficam em `outputs/alula_v040_launch_final/`;
instalação e ensaios exploratórios foram preservados separadamente.

## Limitações e próxima versão

Permanecem necessários estudos de refinamento viscoso, transição/trips, esteiras
e base finita; comparação independente multielemento viscosa e experimental;
condicionamento/robustez em slots estreitos; sensibilidade a versão Node; e
faixas de geometria/alpha diferentes dos casos aqui exercitados. A biblioteca
BL usa curvas/splines e discretizações internas próprias: preservar os nós de
entrada não transforma seus modelos em equivalentes ao HS ou a uma CFD viscosa
geral. Separação forte, confluência de esteiras e estol não foram validados.
As limitações de Joukowski, fuga finita e benchmark Williams anteriores não
foram removidas pelo novo adaptador.

A execução Windows pede cwd curto (limite explícito de 240 unidades UTF-16),
e UNC não foi validado. Ncrit/trips são globais por superfície, não expostos
por elemento nesta versão. A API neutra continua de dois corpos; single-element
é suportado pelo bridge de verificação. Os avisos não são autorização para
ignorar limites, transferir precisão entre malhas ou usar provisórios em rankings.

**Recomendação para v0.4.1:** a integração fornece uma base para uma etapa
delimitada de verificação/refinamento Alula e testes de robustez, com fontes e
controles fixados. Priorizar malhas >=160 por corpo dentro do orçamento de 700,
sensibilidade de transição/esteira e referências independentes adequadas antes
de ampliar fidelidade ou desempenho. Manter casos de um elemento no bridge
enquanto não houver necessidade justificada de alterar o contrato neutro.
Não qualificar otimização aerodinâmica ou validação física nesta entrega.

## Reprodução

```powershell
.\.venv\Scripts\python.exe scripts/install_alula.py --download --portable-node
.\.venv\Scripts\python.exe examples/demo_alula.py --output outputs/alula_reproducao
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe -m pytest -q --reference-data=outputs/reference_data_v011 --williams-data=outputs/reference_data_v032
```

Escolha uma pasta nova/vazia dentro do workspace. Se a execução começar de
outro cwd, passe `--node` e `--alula` absolutos. O relatório não depende de uma
atualização automática do upstream ou de dados baixados durante o solve.
