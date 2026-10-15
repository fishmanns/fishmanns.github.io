# MCP(Model Context Protocol)

## 一瞥

模型上下文协议,简称 **MCP**,是一种开放协议,让 AI 以统一的方式发现和使用外部工具,读取数据以及获取其他上下文信息

如:

+ 读取本地代码文件
+ 查询 SQL 数据库
+ 调用网页或第三方 API
+ 操作 Git 仓库

MCP 的目标是让工具提供者实现标准 **MCP Server**,使得支持 MCP 的 Agent 通过 **MCP CLient** 调用

## 架构

**Host/Client/Server**

**Host**: 用户端实际使用的 Agent 应用

**Client**: 协议客户端,是 Host 内部某个和 **Server** 通信的组件,负责和 Server 发送请求并接收结果,一个 Host 能拥有多个 Client 

**Server**: 能力提供者,实现 MCP 的工具,向 Client 暴露**Tools/Resources/Prompts**等,Server 是**协议角色**,不是物理部署概念,本地也能称为 Server

## 能力

**Tools/Resources/Prompts**

**Tools**: Server 暴露的可执行能力,如搜索文件/查询数据库/调用 API

常见协议方法:

+ tools/list: 查询可用工具
+ tools/call: 调用指定工具

**Resources**: Server 暴露的数据,如配置文件/文档/项目说明等

常见协议方法:

+ resources/list: 查询资源
+ resources/read: 读取指定资源
+ resources/templates/list: 查询资源模板

**Prompts**: Server 暴露的提示模板,用于提供可复用的提示模板或工作流入口,如代码审查模板/日志分析模板

常见协议方法:

+ prompts/list: 查询模板
+ prompts/get: 获取指定模板参数及内容

MCP 并不强制某种工作的实现方式,读取文件既可以设计成 Resources,也可以归类为 Tools

## 协议

MCP 协议基于 JSON-RPC 2.0,一端通过消息告诉另一端要调用什么方法,另一端返回结果,消息媒介通过 json 传播