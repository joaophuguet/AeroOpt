# Instruções permanentes — AeroOpt

- Consulte `docs/PROJECT_CONTEXT.md` antes de alterar a arquitetura.
- Preserve e consulte `docs/GEOMETRY_CONVENTIONS.md` para alterações geométricas.
- Use exclusivamente `.venv\Scripts\python.exe` para operações Python do projeto.
  Não recrie a `.venv` nem instale dependências globalmente.
- Execute os testes após modificações: `.\.venv\Scripts\python.exe -m pytest -q`.
- Atualize a documentação quando decisões técnicas ou a API forem alteradas.
- Não confunda resultados geométricos com resultados aerodinâmicos.
- Não modifique dados de entrada silenciosamente. Não suavize, interpole, reordene
  ou redistribua pontos automaticamente.
- Não introduza funcionalidades fora do escopo sem justificativa documentada.
- Não execute comandos destrutivos sem autorização. Preserve o histórico Git,
  não altere remotos e não crie commits automaticamente.
- Mantenha arquivos gerados, caches e temporários dentro do workspace.
- A v0.1 é exclusivamente geométrica. CFD, solvers aerodinâmicos, NSGA-II e
  modelos substitutos pertencem a versões futuras.
- A v0.1.1 é uma etapa de validação geométrica. Consulte também
  `docs/VALIDATION_V011.md` antes de alterar referências ou tolerâncias.
- Preserve os downloads originais com SHA-256 e proveniência. Não versione
  coordenadas externas sem condições de redistribuição estabelecidas; os
  originais UIUC locais e seus derivados permanecem em `outputs/`.
- Para incluir os testes publicados: `.\.venv\Scripts\python.exe -m pytest -q --reference-data=outputs/reference_data_v011`.
  Sem os dados, as verificações analíticas permanecem executáveis offline;
  nunca invente medições externas.
- Não confunda origem nominal NACA, mínimo x e referência por maior distância
  à fuga. Use índices explícitos quando houver referência conhecida.
- Contato/folga são globais; duplicação e degeneração são verificações
  intrínsecas na escala própria do elemento. Não relaxe tolerâncias para aprovar
  testes sem reproduzir e justificar o defeito.
- A v0.2.0 introduz `aeroopt.solvers`, uma infraestrutura de execução sem solver
  físico. Consulte `docs/SOLVER_INTERFACE.md` e `docs/VALIDATION_V020.md` antes
  de alterar contratos, identidade, estados ou controle de processos.
- Mantenha coeficientes ausentes como `None`/`null`. Adaptadores simulados não
  podem publicar `cl`, `cd` ou `cm`. Saída zero não estabelece convergência.
- Use o runner e `AttemptContext.run_process` para processos externos: argumentos
  sem shell, logs isolados, timeout e encerramento dos descendentes. Preserve
  dados/manifestos e não altere geometrias para satisfazer um adaptador.
- Preserve a separação entre ID científico determinístico e UUID da tentativa.
  Não inclua horários, caminhos de execução ou PIDs no hash da análise.
- O runner exige `min_clearance` explícito. Referência de momento e parâmetros
  de transição/turbulência não recebem valores físicos implícitos.
- Mudanças futuras devem manter os testes geométricos e os novos testes de
  contratos, falhas, caminhos e processos. Não confunda infraestrutura validada
  com aerodinâmica validada.
- A v0.3.0 inclui `PanelSolver`, potencial 2D incompressível invíscido. Consulte
  `docs/PANEL_METHOD.md` e `docs/VALIDATION_V030.md` antes de alterar método,
  sinais, Kutta ou aceitação numérica. Convergência discreta não certifica
  precisão de discretização nem validação experimental. Preserve `cd=None`.
- Consulte também `docs/VALIDATION_V031.md` para a revisão científica. Preserve
  a flag/diagnóstico de fuga finita não verificada; não use esses coeficientes
  para classificar projetos. Referências de pesquisa não são validação viscosa
  nem substituem estudos de refinamento e condicionamento.
- Consulte `docs/VALIDATION_V032.md` e `docs/NUMERICAL_RELIABILITY.md` antes de
  alterar referências independentes ou critérios de ordenação. A_like/B_like
  são casos conformes declarados, não reproduções certificadas dos aerofólios
  de Williams; o critério do benchmark publicado permanece pendente.
- `aeroopt.verification` recebe histórias e limitações explícitas, sem executar
  solvers ou publicar limites rigorosos de erro. Não converter `None`/falhas em
  zero, não qualificar fugas finitas e não interpretar CM maior como projeto
  melhor. O adaptador permanece 0.3.1 enquanto a produção não mudar.

