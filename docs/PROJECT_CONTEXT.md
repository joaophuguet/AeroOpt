# AeroOpt — contexto persistente do projeto

## Propósito e aplicação

O AeroOpt pretende ser uma plataforma modular para otimização aerodinâmica
multifidelidade. A aplicação inicial é uma aeronave cargueira da SAE Aero Design
Regular, tendo a competição East 2027 como alvo. Esse enquadramento é uma
diretriz do projeto, não uma declaração de conformidade com um regulamento
verificado. Nenhum requisito numérico de competição é incorporado na v0.1.

O objetivo científico de longo prazo é maximizar capacidade de carga útil,
sujeita a restrições de decolagem, voo em circuito e estabilidade. A geometria
seccional deve posteriormente se conectar ao desempenho da aeronave completa.

## Configuração e variáveis

O sistema inicial é 2D: elemento principal e flap simples fendido, com perfis
importados independentemente. O flap permanece fixo em todas as fases do voo;
as configurações aqui geradas são candidatos de projeto, não uma programação
de atuação durante o voo. Não há mecanismo, dobradiça móvel ou cálculo 3D.

O vetor de projeto inicial é `[delta_f_deg, dx_f/c_main, dy_f/c_main]`.
A razão `c_f/c_main` é configurável e fixa para uma base de candidatos, não
uma variável inicial de otimização. A posição nominal é a referência global
da borda de ataque do flap antes dos deslocamentos. As convenções detalhadas
estão em `GEOMETRY_CONVENTIONS.md`; elas são parte do contrato científico.

## Decisões arquiteturais da v0.1

- Pacote Python em `src/aeroopt`, instalado em modo editável na `.venv` local.
- NumPy para coordenadas e transformações, Shapely/GEOS para polígonos e
  distâncias, Matplotlib para figuras, pytest para verificação.
- `Airfoil`: originais imutáveis, referências explícitas, normalização uniforme,
  área fechada e proveniência. O contorno é uma poligonal, sem ajuste de curva.
- `FlapParameters`: razão, posição nominal, ângulo em graus e deslocamentos
  absolutos em cordas principais. `MultiElementAirfoil.configure` retorna uma
  nova configuração, calculada sempre a partir das mesmas referências.
- `validate_elements` recebe os contornos globais e produz razões estruturadas.
  Não há reparo automático de geometrias nem mistura com avaliação de solver.
- Exportações separadas por elemento, com precisão suficiente para round trip,
  metadados versionados, parâmetros, unidades, referências e relatório.
- Funções de plotagem retornam objetos Matplotlib e não abrem janelas.
- Sem commits automáticos e sem instalação global. Dados de entrada não são
  alterados. Caches, testes e resultados permanecem no workspace.

## Objetivos e estratégia multifidelidade futura

Métricas futuras: sustentação seccional máxima `cl_max`, arrasto `cd` e momento
`cm` nas condições operacionais; sensibilidade a Reynolds, ângulo de ataque e
transição; robustez a pequenas perturbações geométricas.

A estratégia prevista é iniciar otimização multiobjetivo com NSGA-II e avaliar
MSES como candidato a solver rápido, sujeito a validação. Uma fidelidade maior
deverá usar RANS com modelagem de transição. A classificação de candidatos deve
ser comparada entre fidelidades antes de confiar nela. Modelos substitutos e
alocação adaptativa de fidelidade podem ser incorporados quando houver dados
confiáveis e critérios de erro. A integração com desempenho de aeronave completa
deve converter resultados seccionais em objetivos de carga útil de forma
verificada e com hipóteses explícitas.

Não está decidido que MSES é adequado a toda a região de projeto. Solvers
viscosos rápidos e RANS dependem de modelos, condições de transição, malhas e
convergência. Separação, escoamentos próximos de estol e interação multielemento
exigem validação; nenhuma classificação ou coeficiente deve ser fabricado ou
inferido a partir da aparência geométrica. Resultados 2D não representam
diretamente asa finita, estabilidade ou desempenho da aeronave.

## Escopo implementado e verificação

A v0.1 cobre DAT estrito, normalização de contornos TE-LE-TE, referências
amostradas ou LE explícita, escala/rotação/translação, validação planar,
distância mínima global, exportação JSON/DAT e figuras PNG/SVG. Testes usam
diamantes e quadrados com distâncias/áreas conhecidas, além de duas poligonais
de demonstração. O exemplo registra aceitação e rejeição de seis candidatos.

Nenhum solver aerodinâmico, CFD, NSGA-II ou modelo substituto está implementado.
A execução e seus resultados são registrados em `docs/VALIDATION.md` após os
testes e o exemplo. Os intervalos de dependências e o lock do ambiente validado
permitem separar a API suportada do ambiente efetivamente exercitado.

## Próximas etapas e decisões pendentes

1. Ampliar a validação com perfis reais e condições de integração: a v0.1.1
   verificou NACA 0012, NACA 4412 e S1223, mas não garante identificação automática
   correta da referência nominal nem qualidade suficiente para qualquer solver.
2. Definir gap/overlap aerodinâmicos por referências, eixos e regiões de busca
   explícitos; validar contra convenções do solver e interpretação experimental.
3. Adicionar catálogos/proveniência com hashes de entrada, esquemas de metadados
   e reconstrução de montagens a partir de JSON; tornar exportação transacional.
4. Avaliar discretização controlada e verificações de qualidade sem alterar
   silenciosamente os contornos fornecidos.
5. Definir interface de solver, condições de operação, rastreamento de convergência
   e conjuntos de referência. Somente então iniciar comparação de fidelidades.
6. Validar regulamento aplicável, modelos de voo e restrições do problema antes
   de formular a otimização da aeronave.

Também estão pendentes política para perfis com LE não amostrada, referências
TE fora das extremidades do arquivo, incerteza de discretização e tratamento de
arquivos incompatíveis via conversores explícitos. A posição nominal usada no
exemplo não é uma recomendação aerodinâmica.

## Incremento v0.1.1 — validação geométrica

Foram exercitados três arquivos UIUC preservados e identificados por SHA-256,
além de NACA 0012/4412 gerados analiticamente. A permissão de redistribuição das
coordenadas não foi estabelecida: originais e derivados ficam em `outputs/`,
ignorados pelo Git. Apenas o catálogo de fontes/hashes faz parte do código.
O exemplo `examples/validate_real_airfoils.py` funciona offline com referências
analíticas ou, explicitamente, com downloads/arquivos locais fixados por hash.

Uma implementação escalar independente das equações fornece referências contínuas
por integração e pelo mínimo da curva de nariz em um caso de tampa de fuga.
O estudo distingue erro poligonal, referência de normalização e erro de
implementação. Duas falhas da v0.1 foram reproduzidas: área absoluta rejeitava
flaps pequenos resolvidos, e proximidade absoluta de pontos fazia uma subdivisão
legítima de segmento alterar a validade. Na v0.1.1, qualidade intrínseca do
elemento usa sua própria escala, enquanto contato/folga mantêm a corda principal.
As constantes numéricas não foram relaxadas.

A referência automática por ponto mais distante permanece compatível, mas não
coincide com a origem nominal do NACA 4412 em discretizações mais finas. Casos
com referência conhecida usam índice explícito. Acordo entre duas discretizações
não é certificado de exatidão: há um platô mensurável de erro de distância entre
161 e 321 pontos no estudo. Resultados completos estão em `VALIDATION_V011.md`.

Este incremento continua exclusivamente geométrico. Há evidência suficiente
para iniciar um adaptador de solver sob contratos e referências revisados;
não há validação de painelamento, convergência ou classificação aerodinâmica.
CFD, MSES, RANS, XFOIL e otimização permanecem pendentes.

## Incremento v0.2.0 — interface e infraestrutura de execução

`aeroopt.solvers` separa condições operacionais, contrato abstrato, snapshot de
entrada, execução, interpretação e persistência. Não contém solver físico.
O adaptador simulado é explicitamente não físico e mantém todos os coeficientes
ausentes. Versão instalada: 0.2.0, sem novas dependências de runtime.

O ID SHA-256 usa coordenadas globais, condições, identidade/configurações do
adaptador, esquema e folga explícita. Tentativas têm UUID e arquivos próprios;
pastas/timestamps não entram na identidade científica. Float64 é codificado
exatamente no hash, sem arredondamento de coordenadas. Objetos e snapshots são
imutáveis. O resultado preserva proveniência, condições, estados separados de
execução/convergência, logs relativos e hashes dos artefatos.

Reynolds é baseado na corda principal. Alpha é expresso em graus e define
`V_inf/|V_inf|=(cos(alpha),sin(alpha))` no quadro global; o runner não gira a
geometria. Adaptadores devem converter convenções explicitamente. `cm` presente
exige posição, sentido e normalização declarados, sem referência implícita.

O runner valida geometria/folga antes de preparar o adaptador, recusa alterações
nos dados globais e captura falhas por etapa. Processos externos usam argumentos
sem shell, saída em arquivos e contenção por Job Object no Windows com um worker
bloqueado até a associação. Timeout e saída normal encerram descendentes comuns.
Chamadas Python em processo são cooperativas, sem preempção arbitrária. A API
de caminhos contém os artefatos, mas não constitui um sandbox para código hostil.

Essa infraestrutura permite iniciar um adaptador físico posteriormente, com
faixas suportadas, modelos de transição, critérios de convergência e convenções
validados de forma independente. Não se implementou MSES, XFOIL, RANS,
otimização ou modelos substitutos. Validade geométrica e sucesso do processo
não demonstram fidelidade física. Detalhes e evidências estão em
`SOLVER_INTERFACE.md` e `VALIDATION_V020.md`.

